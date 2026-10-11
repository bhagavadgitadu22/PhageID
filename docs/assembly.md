# Assembly and selection

The pipeline selects one assembly per sample using the highest CheckV quality among its retained viral contigs. Rescue is conditional: it runs only when an earlier attempt does not produce a `Complete` or `High-quality` viral contig. All retained viral contigs from the selected assembly are kept, including those with lower quality.

Paths below are relative to `results_dir`, except workflow source files. Run Snakemake commands from `workflow/`.

## Reads and host filtering

The input `data/samples.tsv` has no header: sample name, long-read input, single-end TruSeq input, and optional bacterial host genome. A read input can be a FASTQ file or a directory of FASTQ files, including gzip-compressed files. If both read types are supplied, the long-read route is used; this is not a hybrid assembly.

In [`01_preprocessing_reads.smk`](../workflow/rules/01_preprocessing_reads.smk):

| Rule | Role |
| --- | --- |
| `concat_reads` | Combine files for the selected read type. Empty input is an error. |
| `preprocess_reads_porechop` | Trim adapters from long reads using Porechop. |
| `preprocess_reads_fastp` | Process single-end short reads with fastp, using quality threshold 20 and minimum length 50 bp. |
| `remove_bacterial_contamination` | Map reads to the supplied host and retain unmapped reads. Without a host genome, copy all preprocessed reads. |
| `assembly_read_status` checkpoint | Use [`assembly_read_stats.py`](../workflow/scripts/assembly_read_stats.py) to count preprocessed and host-filtered reads and record the percentage removed. |

The read sets used by assembly are:

- `filtered`: bacterial-free reads, used first when a host genome is supplied and at least one read remains.
- `all`: all **preprocessed** reads, including host-mapping reads. This does not mean the original untrimmed FASTQs.

If no host genome is supplied, or filtering removes every read, assembly starts directly with `all`. All-read rescue can recover sequences that matched the bacterial reference, including possible prophages.

## Primary assembly

Rules are in [`02_assembly.smk`](../workflow/rules/02_assembly.smk).

### Long reads: Flye

`subsample_reads_flye` runs [`subsample_flye_reads.py`](../workflow/scripts/subsample_flye_reads.py) before `assembly_reads_flye`.

The rule parameters assume a 100,000 bp genome and aim for 1,000× coverage: approximately **100 Mb of read bases**. Smaller read sets are retained. Sampling uses seed 42 and writes a subsampling report; a whole final read can take the total slightly above the target. This is a sampling budget, not an estimate of the actual phage genome length.

Flye uses `--meta --nano-raw`. The rule requires its assembly FASTA, assembly graph, and assembly information table. Full preprocessed reads remain available for mapping and the alternative assembly route.

### Short reads: SPAdes

If no long reads are supplied, `assembly_reads_spades` uses SPAdes with `-s` for single-end input and `--isolate`. It copies `contigs.fasta` to the common candidate `assembly.fasta` path. Autocycler is not used for short-read-only samples.

Both primary assemblers run through [`run_assembly_attempt.py`](../workflow/scripts/run_assembly_attempt.py), which records success or failure without relying on tool-specific error messages.

## Candidate assessment

Every attempted assembly is stored separately:

```text
<sample>/assembly_attempts/<filtered|all>/<flye|spades|autocycler>/
```

The rules in [`04_assembly_selection.smk`](../workflow/rules/04_assembly_selection.smk) assess each candidate before deciding whether rescue is needed:

| Rule | Role |
| --- | --- |
| `genomad_candidate` | Detect viral sequences and provide taxonomy/provirus metadata. An empty failed assembly skips geNomad and produces empty outputs with the necessary table header. |
| `keep_long_viral_contigs_candidate` | Retain viral sequences at least 2,000 bp long. |
| `rename_contigs_candidate` | Prefix contig IDs with the sample name to avoid collisions between samples. |
| `fix_circular_viral_contigs_per_sample_candidate` | Run [`correct_phage_contigs.py`](../workflow/scripts/correct_phage_contigs.py) to correct likely concatemers and terminal repeats, producing corrected FASTA and two CSV reports. |
| `checkv_candidate` | Assess the corrected viral sequences. Recreate the candidate's CheckV directory on each scheduled run so old gene predictions cannot be reused with a changed FASTA. Empty viral input skips CheckV and writes a header-only quality table. |
| `candidate_viral_quality` checkpoint | Expose the completed quality table and assembly status to selection. |

Correction therefore happens **before CheckV and selection**. No viral contigs is a valid assessment outcome; a CheckV, geNomad, or correction execution error still fails its own rule.

## Rescue and final choice

Selection functions are defined in [`Snakefile`](../workflow/Snakefile): `candidate_choice`, `selected_candidate_inputs`, and `selection_manifest_inputs`.

An attempt is ranked by its best retained viral contig:

```text
Complete > High-quality > Medium-quality > Low-quality > Not-determined
         > successful assembly with no retained viral contigs
         > failed assembly attempt
```

There is no combined numerical quality score, and contigs are not pooled across assemblies.

### Long-read route

1. Run Flye on the initial read set and assess it.
2. If it has any `Complete` or `High-quality` viral contig, keep Flye and do not run Autocycler for that read set.
3. Otherwise run Autocycler on the same read set and assess its consensus. Keep whichever has the higher best-contig quality. Ties keep Flye.
4. If the initial read set was `filtered` and its best result is below `High-quality`, repeat this Flye/Autocycler procedure with `all` reads.
5. Compare the best filtered and all-read results. Select the all-read result only if it is strictly better; ties keep the filtered result.

For example, filtered Flye with `Low-quality`, filtered Autocycler with `Medium-quality`, and all-read Flye with `High-quality` selects all-read Flye. Conversely, filtered Flye with `High-quality` stops rescue even if an untested alternative might contain a complete contig.

### Short-read-only route

1. Run SPAdes on the initial read set and assess it.
2. If that set is `filtered` and no `Complete` or `High-quality` viral contig is found, retry SPAdes with `all` reads.
3. Keep the result with the higher best-contig quality. Ties keep the filtered result.

Ordinary assembly failures are ranked below completed assemblies and allow these rescue steps to proceed. Resource and configuration errors are handled separately, as described below.

## Autocycler internals

Rules are in [`03_autocycler_rescue_cycle.smk`](../workflow/rules/03_autocycler_rescue_cycle.smk).

| Rule or script | Role |
| --- | --- |
| `autocycler_genome_size` | Use `autocycler_genome_size` from configuration. The default `auto` invokes the helper's estimation; a failed estimation falls back to 100,000 bp. |
| `autocycler_subsample` | Generate subsets from the current read set. `params.subset_count` controls their number; it is currently four. Ordinary subsampling failure creates a rejection marker and clears partial subset FASTQs. |
| `autocycler_candidate_jobs` | Expand the assembler × subset combinations using the same subset-count parameter. |
| `autocycler_candidate_assembly` | Run one helper assembler on one subset in a separate job, through `run_assembly_attempt.py`. A subsampling rejection skips the helper and records a failed attempt. |
| `autocycler_candidate_assemblies` | Use [`collect_autocycler_candidates.py`](../workflow/scripts/collect_autocycler_candidates.py) to collect successful candidate FASTAs and write `candidate_status.tsv`. |
| `assembly_reads_autocycler` | Run [`autocycler_consensus.sh`](../workflow/scripts/autocycler_consensus.sh) through the assembly wrapper. |

The eight helper assemblers are Canu, Flye, metaMDBG, miniasm/Minipolish, NECAT, NextDenovo/NextPolish, Plassembler, and Raven. With four subsets this gives 32 candidate jobs per attempted read set. Subsets are chosen by Autocycler; four subsets does not imply exactly 25% of the reads in each.

The helper read type is `ont_r10`. Plassembler additionally receives `--args --no_chromosome`, with `PLASSEMBLER_DB` exported from the rule parameter. Other helpers use their standard arguments.

Consensus requires at least two successful candidate assemblies. The script runs compress, cluster, trim and resolve for passing clusters, then combine. Too few successful candidates, no clusters passing QC, or another ordinary command failure records a failed consensus and permits the selection policy to continue.

There is no extra pipeline minimum-read threshold. For example, Autocycler's own “input reads are too shallow to subset” rejection becomes a failed filtered attempt, after which selection can try all reads.

## Ordinary failures versus execution errors

An ordinary failure is **not proof of a biological problem**. The pipeline records `assembler_failed` because the same exit code could also come from a software or parameter issue; its log remains available for investigation.

| Event | Handling |
| --- | --- |
| Assembler returns an ordinary nonzero exit code | Discard partial assembly outputs, record `failed`, and allow rescue. No error-message parsing. |
| Assembler exits zero but required outputs are missing or empty | Record a failed attempt and allow rescue. |
| Flye cannot find long enough reads or assemble disjointigs | Ordinary assembler failure; allow rescue. |
| Autocycler subsampling fails ordinarily, or consensus has no usable clusters | Ordinary failed attempt; allow rescue. |
| Executable cannot be launched, or shell reports exit 126/127 | Execution/configuration error; fail the job. |
| Command dies from a signal, including a shell's signal-style exit code | Execution error; fail the job. |
| Slurm reports OOM, timeout, node failure, cancellation, preemption, or another failed job/step | Execution error; fail the job, even if the wrapper returned successfully. |
| Assessment, source validation, or publication fails | Fail the rule rather than publish incomplete or inconsistent selected results. |

### Enable Slurm checks

The scheduler distinction requires [`slurm_cluster_status.py`](../workflow/scripts/slurm_cluster_status.py). Add this to the existing Snakemake cluster command, from `workflow/`:

```bash
--cluster-status "python scripts/slurm_cluster_status.py"
```

For example:

```bash
SMK_CLUSTER="sbatch -q {cluster.qos} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time}"

snakemake -j 200 --local-cores 1 --use-conda --conda-frontend conda \
  --keep-going --rerun-incomplete \
  --cluster-config ../config/slurm.yaml --cluster "$SMK_CLUSTER" \
  --cluster-status "python scripts/slurm_cluster_status.py" \
  -p viral_report_global
```

The checker queries `sacct` for the parent job **and its steps**. Any failed step overrides a successful parent state, particularly when an assembler subprocess was killed but its wrapper survived. Completed jobs return `success`; active or not-yet-accounted jobs return `running`. Accounting query errors are real errors, not assumed success. The controller must have access to `sacct`.

Without this option, Snakemake's usual completion markers cannot reliably identify a scheduler failure hidden by a surviving wrapper. Local runs detect direct command interruptions but have no Slurm accounting. Memory allocations remain the responsibility of the submission setup; this change adds no Flye memory request.

`--keep-going` permits independent samples to progress after a job error. It does not make the global report complete while one of its dependencies is still failing.

## Publication, reports, and checkpoints

`assembly_selection_manifest` is an **ordinary local rule**, not a checkpoint. It waits for selection's upstream checkpoints to resolve, then uses [`assembly_selection_manifest.py`](../workflow/scripts/assembly_selection_manifest.py) to write `assembly_selection/selected_sources.json`. The manifest records the chosen assembler, read set, and exact source paths. It also tracks the quality/status files of assessed alternatives so changes can trigger reselection.

`select_viral_assembly` is another ordinary local rule. It reads the manifest, validates all sources, and computes read statistics before publishing any selected outputs. Its outputs include:

- `assembly_selection/assembly.fasta`: original selected assembly.
- `circularisation/circular_viruses.fasta`: retained, corrected viral contigs.
- Selected CheckV/geNomad tables, correction reports, assembler/status files, and `assembly_selection/read_stats.tsv`.

In [`05_assembly_report.smk`](../workflow/rules/05_assembly_report.smk), the two mapping rules map all and bacterial-free reads to the corrected viral contigs using theBIGbam `--circular`. `mapped_assembly_info` produces lengths, both coverage values, and Flye circularity metadata. `viral_report` uses [`viral_contig_report.py`](../workflow/scripts/viral_contig_report.py) to write per-sample contig and summary tables.

The `viral_report_global` checkpoint concatenates those into `reports/viral_contigs/viral_contigs.tsv` and `samples.tsv`. This generated `samples.tsv` is separate from the input sample list. Its sample statuses are:

| Status | Meaning |
| --- | --- |
| `ready` | At least one retained viral contig; annotation and comparison are enabled. High CheckV quality is not required for this status. |
| `no_viral_contigs` | An assembly completed, but no viral contigs survived selection. Keep its summary row and skip downstream annotation. |
| `assembly_failed` | All assessed assembly alternatives failed ordinarily. Keep its summary row and skip downstream annotation. This is not a biological diagnosis. |

An unresolved execution error prevents the final report from completing for that sample; it is not silently converted to `assembly_failed`.

Both reports use underscore-separated headers. `Number_reads` counts all preprocessed reads, `Number_reads_without_host` counts retained host-filtered reads, `Assembler` identifies the selected tool, `Assembly_method` indicates `all` or `bacterial free`, and `Assembly_used_reads` counts the actual assembly input after any Flye subsampling. Without a host reference, the first two counts are equal. Sample metadata remains available even when no viral contigs are retained.

### Attempt diagnostics and retries

The wrapper writes `<attempt-log>.attempt.json` before launching its command, preserving the Slurm job ID, log path, execution status, outcome, and exit code. This permits diagnosis even when a whole job is killed before cleanup.

[`assembly_attempt_report.py`](../workflow/scripts/assembly_attempt_report.py) combines these records and Slurm accounting into `reports/viral_contigs/assembly_attempts.tsv` when the global report completes. It includes individual Autocycler helper attempts. While errors block the global report, generate a separate diagnostic report:

```bash
python scripts/assembly_attempt_report.py \
  --results-dir /scratch/boutroux/AKIRA \
  --output /scratch/boutroux/AKIRA/reports/assembly_attempts.tsv
```

Unavailable accounting remains explicitly unknown. Older attempt outputs without metadata are marked `legacy_unknown`; the script cannot reconstruct their cause from a status file alone.

Completed ordinary failures are valid workflow outputs and are normally reused. Adding samples does not automatically retry old failed attempts. `--rerun-incomplete` retries incomplete execution failures; forcing `select_viral_assembly` only republishes selection. To intentionally retry an old assembly, force its assembler rule while targeting a file for the relevant sample, then let selection update. A selected sample target is `<results_dir>/<sample>/assembly_selection/read_stats.tsv`; the directory itself is not a rule output.
