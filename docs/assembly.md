# Assembly and selection

The pipeline selects one assembly per sample using the highest CheckV quality among its retained viral contigs. Rescue steps are conditional: they run only when an earlier attempt does not produce a `Complete` or `High-quality` viral contig. All retained viral contigs from the selected assembly are kept, including those with lower quality.

Paths below are relative to `results_dir`, except workflow source files. Run Snakemake commands from `workflow/`.

## Reads and host filtering

The input `data/samples.tsv` has no header: sample name, long-read input, single-end TruSeq input, and optional bacterial host genome. Read and host paths may be absolute or relative to `sample_data_root` in `config/config.yaml`. The included study sheet uses relative paths. A read input can be a FASTQ file or a directory of FASTQ files, including gzip-compressed files. If both read types are supplied, the long-read route is used; this is not a hybrid assembly.

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
- `all`: all **preprocessed** reads, including host-mapping reads.

If no host genome is supplied, or filtering removes every read, assembly starts directly with `all`. All-read rescue can recover sequences that matched the bacterial reference, including possible prophages.

## Primary assembly

Rules are in [`02_assembly.smk`](../workflow/rules/02_assembly.smk). Both short and long read primary assemblers run through `run_assembly_attempt.py`. 

### Long reads: Flye

`subsample_reads_flye` runs [`subsample_flye_reads.py`](../workflow/scripts/subsample_flye_reads.py) before `assembly_reads_flye`. The rule parameters assume a 100,000 bp genome and aim for 1,000× coverage: approximately **100 Mb of read bases**. Smaller read sets are retained. Flye uses `--meta --nano-raw`.

### Short reads: SPAdes

If no long reads are supplied, `assembly_reads_spades` uses SPAdes with `-s` for single-end input and `--isolate`. It copies `contigs.fasta` to the common candidate `assembly.fasta` path. Autocycler is not used for short-read-only samples.

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
| `checkv_candidate` | Assess the corrected viral sequences. |
| `candidate_viral_quality` checkpoint | Expose the completed quality table and assembly status to selection. |

## Rescue and final choice

Several assembly steps are tried depending on the read set (short or long) until a good-enough assembly is obtained. Selection functions are defined in [`Snakefile`](../workflow/Snakefile): `candidate_choice`, `selected_candidate_inputs`, and `selection_manifest_inputs`.

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

## Autocycler internals

Rules are in [`03_autocycler_rescue_cycle.smk`](../workflow/rules/03_autocycler_rescue_cycle.smk).

| Rule or script | Role |
| --- | --- |
| `autocycler_genome_size` | Use `autocycler_genome_size` from configuration. The default `auto` invokes the helper's estimation; a failed estimation falls back to 100,000 bp. |
| `autocycler_subsample` | Generate subsets from the current read set. `params.subset_count` controls their number (default 4). |
| `autocycler_candidate_jobs` | Expand the assembler × subset combinations using the same subset-count parameter. The eight helper assemblers are Canu, Flye, metaMDBG, miniasm/Minipolish, NECAT, NextDenovo/NextPolish, Plassembler, and Raven. With four subsets this gives 32 candidate jobs per attempted read set. |
| `autocycler_candidate_assembly` | Run one helper assembler on one subset in a separate job, through `run_assembly_attempt.py`. |
| `autocycler_candidate_assemblies` | Use [`collect_autocycler_candidates.py`](../workflow/scripts/collect_autocycler_candidates.py) to collect successful candidate FASTAs and write `candidate_status.tsv`. |
| `assembly_reads_autocycler` | Run [`autocycler_consensus.sh`](../workflow/scripts/autocycler_consensus.sh) through the assembly wrapper. |

Consensus requires at least two successful candidate assemblies. The script runs compress, cluster, trim and resolve for passing clusters, then combine. Too few successful candidates, no clusters passing QC, or another ordinary command failure records a failed consensus and permits the selection policy to continue.

## Ordinary failures versus execution errors

Assembly errors need to be distinguished from jobs that were interrupted midway. Flye can fail because there are too few long reads to assemble disjointigs, for example, which is an acceptable error. Conversely, an assembly job failing because of a Slurm error (OOM out-of-memory, timeout, node failure, cancellation, preemption, or another failed job/step) means the assembly did not finish its job but could work if given more resources. 

To distinguish ordinary failed jobs from interrupted ones, the Snakemake run requires [`slurm_cluster_status.py`](../workflow/scripts/slurm_cluster_status.py):

```bash
SMK_CLUSTER="sbatch -q {cluster.qos} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time}"

snakemake -j 200 --local-cores 1 --use-conda --conda-frontend conda \
  --keep-going --rerun-incomplete \
  --cluster-config ../config/slurm.yaml --cluster "$SMK_CLUSTER" \
  --cluster-status "python scripts/slurm_cluster_status.py" \
  -p viral_report_global
```

The checker queries `sacct` for the parent job **and its steps**. The wrapper writes `<attempt-log>.attempt.json` before launching its command, preserving the Slurm job ID, log path, execution status, outcome, and exit code. This permits diagnosis even when a whole job is killed before cleanup.

## Final assembly selection

`assembly_selection_manifest` waits for selection's upstream checkpoints to resolve, then uses [`assembly_selection_manifest.py`](../workflow/scripts/assembly_selection_manifest.py) to record the chosen assembler, read set, and exact source paths. It also tracks the quality/status files of assessed alternatives so changes can trigger reselection.

Another rule, `select_viral_assembly`, reads the manifest, validates all sources, and computes read statistics before publishing the selected automatic outputs. Its outputs include:

- `assembly_selection/assembly.fasta`: original selected assembly.
- `assembly_selection/automatic_circular_viruses.fasta`: automatically corrected viral contigs.
- Selected CheckV/geNomad tables, automatic correction reports, assembler/status files, and `assembly_selection/read_stats.tsv`.

You can correct manually this automatic assembly if needed. It can be useful to remove a repeat or make another sequence edit for example. To do that, copy `assembly_selection/automatic_circular_viruses.fasta`, correct it (you need to keep the same contig IDs though) and put the corrected FASTA to the stable manual input path:

```bash
results_dir=$(python -c 'import yaml; print(yaml.safe_load(open("../config/config.yaml"))["results_dir"])')
sample="$results_dir/Brevumdimonas307_T9_22"
mkdir -p "$sample/manual_corrections"
cp "<YOUR_VIRAL_GENOME>.fasta" "$sample/manual_corrections/circular_viruses.fasta"
```

The rule `finalize_viral_assembly` publishes the final files for a sample. If no manual override was done, it will use the outputs from the `select_viral_assembly` rule to publish `circularisation/circular_viruses.fasta`, its correction reports, and `checkv/quality_summary.tsv`. With an override, it copies the manual viral genome to `circularisation/circular_viruses.fasta` and reruns CheckV on the final sequence. The separate `circularisation/manual_corrections_report.tsv` records which sequences changed and their length differences.

Editing the manual input at any time triggers the rule `finalize_viral_assembly` again. To remove an override, delete the file under `manual_corrections/` and rerun Snakemake.

## Final assembly reports

In [`05_assembly_report.smk`](../workflow/rules/05_assembly_report.smk), the per-sample reports and global sample summary mentioned in [the README](../README.md) contain a column Status explaining the assembly status. Three values are possible:
- `ready`: At least one retained viral contig; annotation and comparison are enabled. High CheckV quality is not required for this status.
- `no_viral_contigs`: An assembly completed, but no viral contigs survived selection. Keep its summary row and skip downstream annotation.
- `assembly_failed`: All assessed assembly alternatives failed ordinarily. Keep its summary row and skip downstream annotation. This is not a biological diagnosis.

An unresolved execution error prevents the final report from completing for that sample; it is not silently converted to `assembly_failed` (see [Ordinary failures versus execution errors](#ordinary-failures-versus-execution-errors)).

A report showing the outcome of each assembly attempt is made by [`assembly_attempt_report.py`](../workflow/scripts/assembly_attempt_report.py) into `reports/assembly_attempts.tsv` when the global report completes. It includes individual Autocycler helper attempts. If errors block the global report, you can generate a separate diagnostic report with:

```bash
results_dir=$(python -c 'import yaml; print(yaml.safe_load(open("../config/config.yaml"))["results_dir"])')
python scripts/assembly_attempt_report.py --results-dir "$results_dir" --output "$results_dir/reports/assembly_attempts.tsv"
```
