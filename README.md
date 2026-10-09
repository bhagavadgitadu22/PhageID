# ALP Snakemake workflow for phage genome characterization

PhageID assembles phage isolates from Oxford Nanopore long reads or TruSeq single-end short reads, selects viral contigs, assesses their quality, annotates them and compares samples.

## Assembly and checkpoints

1. **Preprocess reads.** Combine FASTQ inputs, then trim long reads with Porechop or process short reads with fastp. If a bacterial host genome is supplied, remove host-mapping reads from either read type. The `assembly_read_status` checkpoint records read counts and the percentage removed..
2. **Assemble.** Start with bacterial-free reads when available; otherwise use all preprocessed reads. Long-read samples use Flye after subsampling to approximately 1,000× (hypothesizing a 100 kbp assembly). Short-read-only samples use SPAdes (`--isolate`, single-end input). When both read columns are provided, the long-read route is used.
3. **Assess each candidate.** Run geNomad, retain viral sequences of at least 2 kbp, prefix contig IDs with the sample name, correct likely concatemers and terminal repeats, then run CheckV. The `candidate_viral_quality` checkpoint makes these results available for assembly selection.
4. **Rescue with Autocycler.** If Flye has no `Complete` or `High-quality` viral contig, run Autocycler on the same read set. Four subsets are assembled independently with eight assemblers: Canu, Flye, metaMDBG, miniasm/Minipolish, NECAT, NextDenovo/NextPolish, Plassembler and Raven, followed by consensus assembly. Compare the highest individual CheckV quality found in Flye and Autocycler.
5. **Rescue with all reads and select.** If neither filtered assembly has a `Complete` or `High-quality` viral contig, repeat steps 2–4 with all preprocessed reads, including bacterial reads, to allow host-matching prophage recovery. Short-read-only samples likewise retry SPAdes with all reads when the filtered assembly has no `Complete` or `High-quality` viral contig. Select the attempted assembly with the highest individual CheckV category: `Complete` > `High-quality` > `Medium-quality` > `Low-quality` > `Not-determined` > no viral contigs. Keep **all** retained corrected viral contigs from that assembly, including lower-quality ones; contigs from different assemblies are not merged. Ties prefer the earlier attempt: filtered reads before all reads, and Flye before Autocycler within a read set. Once a read set yields a `Complete` or `High-quality` contig, further rescue is unnecessary. An assembler command failure still stops the workflow; rescue follows completed quality assessments.
6. **Report and gate downstream work.** Map all reads and bacterial-free reads onto the selected corrected viral contigs using theBIGbam `--circular`. Write per-sample contig and summary tables. The `viral_report_global` checkpoint concatenates them and allows only samples with retained viral contigs to enter annotation and comparison. Empty viral results remain in the sample report and skip CheckV execution and downstream annotation.

Autocycler uses `autocycler_subsample.params.subset_count` (currently 4) for both subsampling and candidate jobs. `autocycler_genome_size` in `config/config.yaml` defaults to `auto`, estimated with Raven; a fixed value such as `100000` is also accepted. Candidate read chemistry is currently `ont_r10` in `workflow/rules/03_autocycler_rescue_cycle.smk`. Flye subsampling parameters are set in `workflow/rules/02_assembly.smk`. Completed candidate failures are logged and skipped; consensus requires at least two successful assemblies and one QC-pass cluster. Candidate status and consensus intermediates are retained under each assembly attempt.

Plassembler candidates run in chromosome-free mode (`--no_chromosome`, forwarded through Autocycler `--args`) and use `/work/river/Databases/plassembler_db`, exported as `PLASSEMBLER_DB` from the candidate rule's `params.plassembler_db`. Download the database there, or edit that parameter; the workflow does not download it. See the [Autocycler helper implementation](https://github.com/rrwick/Autocycler/blob/v0.7.0/src/helper.rs) for database lookup and candidate dependencies.

## After assembly

Selected contigs receive Pharokka/PHANOTATE and Phold annotations, alternative Genotate gene calls with Pharokka protein annotation, per-contig plots, PhageTerm analysis, antiDefenseFinder predictions, Empathi/Sublyme predictions and CheckAMG analysis. Functional predictions are added to the Phold GFF for theBIGbam.

The separate comparison target calculates all-versus-all ANI and dereplicates contigs at ≥95% ANI and ≥85% aligned fraction of the member contig, with no minimum representative aligned fraction. Clustering considers contigs in descending length, so representatives are chosen by length. It produces a dereplication report, LoVis4u comparisons and vConTACT3 clustering with prokaryotic reference genomes.

The first theBIGbam database uses MAG view: one sample contains all its selected viral contigs, with two BAMs per sample (all reads and bacterial-free reads). The `prepare_post_dereplication_mapping` checkpoint resolves sample–representative pairs after clustering. The second database uses contig view, representative-only annotations and `{sample}_on_{representative}.bam` mappings.

## Inputs and running

`data/samples.tsv` is tab-separated **without a header**, with columns in this order:

| Sample | Long reads | TruSeq single-end reads | Bacterial host genome |
| --- | --- | --- | --- |
| isolate_long | /path/to/long_reads/ | NA | /path/to/host.fasta |
| isolate_short | NA | /path/to/single_end.fastq.gz | NA |

Each read input can be a file or a directory. Directories combine `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz` files in filename order, ignoring subdirectories. Empty fields, `NA`, `None` and `-` mean absent. At least one read input is required. An absent host genome keeps all reads; a supplied host FASTA must exist.

Set `samples_file` and `results_dir` in `config/config.yaml`, and check the database paths in the rules. Run from `workflow/` using the Snakemake environment defined in `envs/requirements.yaml`:

```bash
cd workflow
snakemake --use-conda --cores 8 corrected_assembly_all_phages  # assembly and reports
snakemake --use-conda --cores 8 pipeline_all_phages           # add sample annotation
snakemake --use-conda --cores 8 pipeline_comparison_phages    # add cross-sample comparisons
```

`pipeline_all_phages` is the default target. Autocycler can also be requested directly by targeting `<results_dir>/<sample>/assembly_attempts/{filtered,all}/autocycler/assembly.fasta` for a long-read sample.

## Main reports

- `<sample>/assembly_stats.tsv`: one row per selected viral contig, with assembler, assembly read set, preprocessing read counts, host-mapping percentage, original and corrected lengths, both coverage values, Flye circularity, correction counts, CheckV metrics and geNomad taxonomy/provirus status.
- `<sample>/viral_sample_report.tsv`: one sample summary row. Global tables are `reports/viral_contigs/viral_contigs.tsv` and `reports/viral_contigs/samples.tsv`.
- `combined_viruses/dereplication/dereplication_report.tsv`: every contig before dereplication, its selection status, representative, lengths, cluster size, ANI and aligned fractions.

`Total bp` is the viral sequence length before correction; for a provirus it is the extracted region length. `Total corrected bp` and coverage refer to corrected sequences. Coverage counts aligned bases, including origin-spanning alignments, divided by corrected length. `Coverage` uses all preprocessed reads; `Coverage without bacteria` uses host-filtered reads. Without a host genome, both use the same reads and bacterial percentage is `NA`. Flye circularity is available for whole Flye contigs; extracted proviruses and other assemblers report `NA`.

Read counts are measured after preprocessing. `Reads used for assembly number` is the actual Flye subset size, or the input read set size for SPAdes/Autocycler; Autocycler assembles overlapping subsets of that set. CheckV metrics describe corrected contigs, while geNomad annotations describe the original viral sequences.

## Key resources table

The tables follow the Cell Press STAR Methods convention. Software identifiers below reflect the environment specifications and installation rules in this workflow. Dependencies installed indirectly without an explicit version constraint are marked as unpinned.

### Software and algorithms

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
| --- | --- | --- |
| Snakemake | [Snakemake](https://snakemake.github.io/) | v7.32.4 |
| Mamba | [mamba-org](https://mamba.readthedocs.io/) | v1.5.6 |
| Python | [Python Software Foundation](https://www.python.org/) | v3.10.13 (workflow); v3.8.10 (Genotate); v3.9 (PhageTerm); v3.11.5 (Empathi and Sublyme); v3.11 (CheckAMG) |
| Porechop | [Wick; Porechop](https://github.com/rrwick/Porechop) | v0.2.4 |
| fastp | [fastp](https://github.com/OpenGene/fastp) | v0.23.4; single-end TruSeq preprocessing |
| Seqtk | [Li; seqtk](https://github.com/lh3/seqtk) | v1.5 |
| Autocycler | [Autocycler](https://github.com/rrwick/Autocycler) | v0.7.0; automatic long-read fallback |
| miniasm and Minipolish | [miniasm](https://github.com/lh3/miniasm); [Minipolish](https://github.com/rrwick/Minipolish) | miniasm v0.3; Minipolish v0.2.1 |
| Racon | [Racon](https://github.com/lbcb-sci/racon) | v1.5.0 |
| Canu | [Canu](https://github.com/marbl/canu) | v2.3 |
| metaMDBG | [metaMDBG](https://github.com/GaetanBenoitDev/metaMDBG) | v1.4 |
| NECAT | [NECAT](https://github.com/xiaochuanle/NECAT) | v0.0.1_update20200803 |
| NextDenovo | [NextDenovo](https://github.com/Nextomics/NextDenovo) | v2.5.2 |
| NextPolish | [NextPolish](https://github.com/Nextomics/NextPolish) | Unpinned; polishing for the NextDenovo candidate |
| Plassembler | [Plassembler](https://github.com/gbouras13/plassembler) | v1.8.5 |
| Raven | [Raven](https://github.com/lbcb-sci/raven) | v1.8.3 (Autocycler candidate assemblies) |
| Flye | [Flye](https://github.com/mikolmogorov/Flye) | v2.9.6 |
| SPAdes | [SPAdes](https://github.com/ablab/spades) | v4.3.0 |
| geNomad | [Camargo et al.; geNomad](https://github.com/apcamargo/genomad) | v1.12.0 |
| CheckV | [CheckV](https://pypi.org/project/checkv/) | v1.1.1 |
| PhageTermVirome | [PhageTermVirome](https://pypi.org/project/phagetermvirome/4.3/) | v4.3 |
| Pharokka | [Bouras et al.; Pharokka](https://github.com/gbouras13/pharokka) | v1.10.1 |
| Phold | [Bouras et al.; Phold](https://github.com/gbouras13/phold) | v1.3.0 |
| Genotate | [Genotate](https://github.com/deprekate/genotate) | v0.15 |
| Empathi | [Boulay; Empathi](https://huggingface.co/AlexandreBoulay/empathi) | Git revision `7e9cccc00b29e11db1c0bfa066167e8026ef981c` |
| Sublyme | [Rousseau Team; Sublyme](https://github.com/Rousseau-Team/sublyme) | v1.2.2 |
| DefenseFinder | [DefenseFinder](https://github.com/mdmparis/defense-finder) | v2.0.0 (`mdmparis-defense-finder`) |
| MacSyFinder | [MacSyFinder](https://macsyfinder.readthedocs.io/) | v2.1.4 |
| HMMER | [HMMER](http://hmmer.org/) | v3.1 (DefenseFinder environment) |
| CheckAMG | [Anantharaman Lab; CheckAMG](https://github.com/AnantharamanLab/CheckAMG) | Git revision `d29aaef` |
| LoVis4u | [Egorov et al.; LoVis4u](https://github.com/art-egorov/lovis4u) | v0.1.5 |
| vConTACT3 | [vConTACT3](https://vcontact3.readthedocs.io/) | v3.1.6 |
| theBIGbam | [theBIGbam](https://github.com/bhagavadgitadu22/theBIGbam) | v0.7.0 |
| Minimap2 | [Li; minimap2](https://github.com/lh3/minimap2) | v2.31 (Autocycler); otherwise unpinned; used for read mapping |
| SAMtools | [SAMtools](https://www.htslib.org/) | Unpinned; used for alignment processing |
| BLAST+ (BLASTn) | [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) | Unpinned; used for self-alignment and genome comparisons |
| Biopython | [Biopython](https://biopython.org/) | v1.80 (Pharokka environment) |
| TensorFlow | [TensorFlow](https://www.tensorflow.org/) | v2.13.1 (Genotate) |
| PyTorch | [PyTorch](https://pytorch.org/) | v2.6 (Phold); v2.7.1, CUDA 12.8 wheels (Empathi); v2.8.0, CPU wheels (CheckAMG) |
| NumPy | [NumPy](https://numpy.org/) | <2.0 in the PhageTerm environment; exact version unpinned |
| Matplotlib | [Matplotlib](https://matplotlib.org/) | Unpinned; used by custom plotting scripts |
| pyCirclize | [pyCirclize](https://github.com/moshi4/pyCirclize) | Unpinned; used for Genotate contig plots |
| GNU Scientific Library | [GNU GSL](https://www.gnu.org/software/gsl/) | v2.7.0 (Pharokka environment) |
| Poetry | [Poetry](https://python-poetry.org/) | v1.8.5 (PhageTerm installation) |
| PhageID workflow and custom scripts | [This repository](https://github.com/bhagavadgitadu22/PhageID) | `workflow/` and `workflow/scripts/`; use the repository commit corresponding to the analysis |

### Deposited data and reference databases

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
| --- | --- | --- |
| Short-read sequencing data | User-supplied single-end reads | Optional third (TruSeq) column in `data/samples.tsv` |
| Long-read sequencing data | User-supplied sample data | Per-sample read files or directories specified in `data/samples.tsv`; accession numbers not recorded in the workflow |
| Host reference genomes | User-supplied reference data | Optional host genome FASTA paths specified in `data/samples.tsv`; accession numbers not recorded in the workflow |
| geNomad database | [geNomad database](https://portal.nersc.gov/genomad/) | `/work/river/Databases/genomad_db/`; release not recorded |
| CheckV database | [CheckV database](https://portal.nersc.gov/CheckV/) | v1.5; `/work/river/Databases/checkv-db-v1.5/` |
| Pharokka database bundle | [Pharokka](https://github.com/gbouras13/pharokka) | `/work/river/Databases/pharokka_db/`; release not recorded |
| Phold database bundle | [Phold](https://github.com/gbouras13/phold) | `/work/river/Databases/phold_db/`; release not recorded |
| Empathi model files | [Empathi](https://huggingface.co/AlexandreBoulay/empathi) | Git LFS files retrieved at revision `7e9cccc00b29e11db1c0bfa066167e8026ef981c` |
| DefenseFinder models | [DefenseFinder models](https://github.com/mdmparis/defense-finder-models) | Downloaded with `defense-finder update`; model release unpinned |
| CheckAMG annotation database | [CheckAMG](https://github.com/AnantharamanLab/CheckAMG) | `CheckAMG_annotate_db_v1.1_20260316`, as named in the configured database path |
| vConTACT3 reference database | [vConTACT3 databases](https://vcontact3.readthedocs.io/en/latest/databases.html) | `/work/river/Databases/vContact3_db`; release not recorded |
| LoVis4u HMM profiles | [LoVis4u](https://github.com/art-egorov/lovis4u) | Downloaded with `lovis4u --get-hmms`; profile release unpinned |

Database identifiers reflect the configured paths and download commands; local database contents have not been verified. DefenseFinder models and LoVis4u HMM profiles depend on the download date. For reproducible analyses, record database release identifiers or checksums and retain the resolved environment package lists alongside the workflow commit.


## Citing the pipeline

If you use this workflow in your research, please cite this paper:
    Wai Hoe Chin, Martin Boutroux, Akira Harding, Davide Demurtas, Florian Baier, Hannes Peter (2026). Viral isolation reveals novel and diverse phages infecting natural stream biofilms. bioRxiv 2026.03.26.713887. https://doi.org/10.64898/2026.03.26.713887
