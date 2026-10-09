# ALP snakemake workflow for phage genome characterization

This Snakemake workflow allows the characterization of phage genomes via long read sequencing.
It preprocesses the long reads using porechop, assembles each sample with flye and detects viral contigs longer than 2 kbp with geNomad.
Viral quality is assessed with CheckV and contigs are annotated via multiple viral tools (pharokka, phold, empathi, sublyme).
All viral contigs dereplicated using MIUVIG standards and final comparisons are produced across samples using ANI, lovis4u and thebigbam.

## Key resources table

The tables follow the Cell Press STAR Methods convention. Software identifiers below reflect the environment specifications and installation rules in this workflow. Dependencies installed indirectly without an explicit version constraint are marked as unpinned.

### Software and algorithms

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
| --- | --- | --- |
| Snakemake | [Snakemake](https://snakemake.github.io/) | v7.32.4 |
| Mamba | [mamba-org](https://mamba.readthedocs.io/) | v1.5.6 |
| Python | [Python Software Foundation](https://www.python.org/) | v3.10.13 (workflow); v3.8.10 (Genotate); v3.9 (PhageTerm); v3.11.5 (Empathi and Sublyme); v3.11 (CheckAMG) |
| Porechop | [Wick; Porechop](https://github.com/rrwick/Porechop) | v0.2.4 |
| Seqtk | [Li; seqtk](https://github.com/lh3/seqtk) | v1.5 |
| Autocycler | [Autocycler](https://github.com/rrwick/Autocycler) | v0.7.0; automatic long-read fallback |
| miniasm and Minipolish | [miniasm](https://github.com/lh3/miniasm); [Minipolish](https://github.com/rrwick/Minipolish) | Unpinned; Autocycler candidate assembly dependencies |
| Racon | [Racon](https://github.com/lbcb-sci/racon) | Unpinned; Autocycler candidate assembly dependency |
| Canu | [Canu](https://github.com/marbl/canu) | Unpinned; Autocycler candidate assembler |
| metaMDBG | [metaMDBG](https://github.com/GaetanBenoitDev/metaMDBG) | Unpinned; Autocycler candidate assembler |
| NECAT | [NECAT](https://github.com/xiaochuanle/NECAT) | Unpinned; Autocycler candidate assembler |
| NextDenovo | [NextDenovo](https://github.com/Nextomics/NextDenovo) | Unpinned; Autocycler candidate assembler |
| Plassembler | [Plassembler](https://github.com/gbouras13/plassembler) | Unpinned; Autocycler candidate assembler |
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
| Minimap2 | [Li; minimap2](https://github.com/lh3/minimap2) | Unpinned; used for read mapping |
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
| PhageID workflow and custom scripts | [This repository](https://github.com/bhagavadgitadu22/PhageID) | `workflow/` and `scripts/`; use the repository commit corresponding to the analysis |

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

If you use this worflow in your research, please cite this paper:
    Wai Hoe Chin, Martin Boutroux, Akira Harding, Davide Demurtas, Florian Baier, Hannes Peter (2026). Viral isolation reveals novel and diverse phages infecting natural stream biofilms. bioRxiv 2026.03.26.713887. https://doi.org/10.64898/2026.03.26.713887

## Sample sheet

`data/samples.tsv` is tab-separated, without a header: sample name, long reads, TruSeq single-end short reads, and optional bacterial host genome FASTA. Either read column accepts a FASTQ file or a directory containing `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` files. Directory inputs are combined in filename order; other files and subdirectories are ignored. Use empty fields, `NA`, `None`, or `-` for missing inputs. At least one read input is required. If long reads are supplied, they are used for assembly and mapping; otherwise the single-end short reads are assembled with SPAdes (`--isolate`) and mapped with minimap2's short-read preset. An absent host genome skips host filtering; a supplied host FASTA must exist.

Example short-read-only row (empty second and fourth fields):
```text
sample_name		/path/to/single_end.fastq.gz
```

Per-sample `assembly_stats.tsv` and the global `reports/viral_contigs/viral_contigs.tsv` include CheckV gene counts, quality, completeness and contamination, and geNomad provirus status and taxonomy. CheckV describes corrected sequences; geNomad describes the original viral sequences. `viral_report` writes each sample’s contig table and a one-row `viral_sample_report.tsv` after CheckV. The `viral_report_global` checkpoint concatenates these into global contig and sample tables, then selects samples with viral contigs for annotation. CheckV is skipped for empty viral FASTAs. Both coverage columns come from the same circular BAMs used by the first theBIGbam database, mapped onto corrected viral contigs: `Coverage` uses all reads and `Coverage without bacteria` uses host-filtered reads. With no host genome, both use the same reads. Provirus coverage is measured on the corrected extracted sequence; its original reported length is the extracted region. `Flye circularity` retains Flye’s Y/N value for whole contigs and is `NA` for extracted proviruses and other assemblers.

The comparison target also writes `combined_viruses/dereplication/dereplication_report.tsv`: one row per corrected viral contig before dereplication, with length, selection status, representative ID and length, and cluster size (including the representative).

The dereplication report also includes `ANI identity (%)`, `Aligned fraction (%)` of the member contig, and `Representative aligned fraction (%)`. These use the representative-to-member ANI row used in clustering: `pid`, `tcov`, and `qcov`, respectively. Representatives compared with themselves are reported as 100% for all three values. Current clustering requires ANI ≥95% and member aligned fraction ≥85%; representative aligned fraction has no minimum (`min_qcov=0`).

vConTACT3 runs on dereplicated representatives. `genemap_for_vcontact3` extracts Pharokka CDS translations from GenBank files into a combined `vcontact3_inputs/proteins.faa`, writes matching `gene2genome.tsv` (`protein_id`, `genome_id`), and provides `genome_lengths.tsv`. Protein IDs are made unique across contigs. The comparison target runs vConTACT3 with prokaryotic reference genomes and exports Cytoscape, profiles, and completeness results to `vcontact3_results/`.

## Autocycler assembly fallback

Long-read samples first run Flye, geNomad, contig correction and CheckV. If no retained viral contig is `Complete` or `High-quality`, a checkpoint automatically runs Autocycler from Porechop reads and assesses its assembly through the same steps. Downstream analyses use the Autocycler contigs when this fallback is triggered. Short-read-only samples use SPAdes without this fallback. You can also request `<results_dir>/<sample>/autocycler/assembly.fasta` directly. The route estimates genome size, creates four read subsets, assembles each with eight assemblers: Canu, Flye, metaMDBG, miniasm/Minipolish, NECAT, NextDenovo, Plassembler and Raven, then compresses, clusters, trims, resolves and combines QC-pass clusters using the [automated Autocycler workflow](https://github.com/rrwick/Autocycler/wiki/Fully-automated-assembly).

From `workflow/`:

```bash
snakemake --use-conda --cores 8 /scratch/boutroux/AKIRA/SAMPLE/autocycler/assembly.fasta
```

Replace the sample and results path. Optional configuration: `autocycler_subsets` (default 4), `autocycler_genome_size` (default `auto`), and `autocycler_read_type` (default `ont_r10`; alternatives `ont_r9`, `pacbio_clr`, `pacbio_hifi`). Set the read type to match your sequencing chemistry.

Failed candidate assemblies are logged and skipped; the run requires multiple successful assemblies and at least one QC-pass cluster. Intermediate graphs and metrics are retained under `autocycler/consensus_work/attempt_*/autocycler_out/`. The selected assembler is recorded in `assembly_selection/assembler.txt` and the `Assembler` column of per-sample and global contig reports. All assemblers reuse the two theBIGbam `--circular` mappings for coverage before and after host filtering. Mean depth counts aligned M/=/X bases, including those wrapping across the origin, divided by the corrected contig length. Secondary, duplicate, unmapped and QC-failed alignments are excluded; supplementary alignments are retained. Only Flye supplies assembler circularity metadata. A Flye command failure still stops the workflow; this fallback applies after successful Flye processing and CheckV assessment.

## First theBIGbam database

`thebigbam_ALP.db` uses `--view mag`: all selected viral contigs from one sample form one MAG named after that sample. Per-sample enriched GFF files and corrected FASTAs are staged under `thebigbam/mag_inputs/`. Each sample contributes two mappings, `{sample}.bam` (all reads) and `{sample}_without_bacteria.bam` (host-filtered reads). With no host genome, both mappings use the same reads. The database after dereplication keeps its existing contig view.
