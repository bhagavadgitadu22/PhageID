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
| Flye | [Flye](https://github.com/mikolmogorov/Flye) | v2.9.6 |
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
| Long-read sequencing data | User-supplied sample data | Per-sample read directories specified in `data/samples.tsv`; accession numbers not recorded in the workflow |
| Host reference genomes | User-supplied reference data | Optional host genome FASTA paths specified in `data/samples.tsv`; accession numbers not recorded in the workflow |
| geNomad database | [geNomad database](https://portal.nersc.gov/genomad/) | `/work/river/Databases/genomad_db/`; release not recorded |
| CheckV database | [CheckV database](https://portal.nersc.gov/CheckV/) | v1.5; `/work/river/Databases/checkv-db-v1.5/` |
| Pharokka database bundle | [Pharokka](https://github.com/gbouras13/pharokka) | `/work/river/Databases/pharokka_db/`; release not recorded |
| Phold database bundle | [Phold](https://github.com/gbouras13/phold) | `/work/river/Databases/phold_db/`; release not recorded |
| Empathi model files | [Empathi](https://huggingface.co/AlexandreBoulay/empathi) | Git LFS files retrieved at revision `7e9cccc00b29e11db1c0bfa066167e8026ef981c` |
| DefenseFinder models | [DefenseFinder models](https://github.com/mdmparis/defense-finder-models) | Downloaded with `defense-finder update`; model release unpinned |
| CheckAMG annotation database | [CheckAMG](https://github.com/AnantharamanLab/CheckAMG) | `CheckAMG_annotate_db_v1.1_20260316`, as named in the configured database path |
| LoVis4u HMM profiles | [LoVis4u](https://github.com/art-egorov/lovis4u) | Downloaded with `lovis4u --get-hmms`; profile release unpinned |

Database identifiers reflect the configured paths and download commands; local database contents have not been verified. DefenseFinder models and LoVis4u HMM profiles depend on the download date. For reproducible analyses, record database release identifiers or checksums and retain the resolved environment package lists alongside the workflow commit.

## Citing the pipeline

If you use this worflow in your research, please cite this paper:
    Wai Hoe Chin, Martin Boutroux, Akira Harding, Davide Demurtas, Florian Baier, Hannes Peter (2026). Viral isolation reveals novel and diverse phages infecting natural stream biofilms. bioRxiv 2026.03.26.713887. https://doi.org/10.64898/2026.03.26.713887

## Sample sheet

`data/samples.tsv` is tab-separated, without a header: sample name, reads directory, and an optional host genome FASTA path. Omit the third field or leave it empty (also accepted: `NA`, `None`, or `-`) when no host genome is available. In that case, host-read filtering is skipped and all adapter-trimmed reads are retained for theBIGbam. Assembly always uses adapter-trimmed reads. A supplied host FASTA path must exist.

Per-sample `assembly_stats.tsv` and the global `reports/viral_contigs/viral_contigs.tsv` include CheckV gene counts, quality, completeness and completeness method, and geNomad provirus status and taxonomy. CheckV describes corrected sequences; geNomad describes the original viral sequences. A single reporting checkpoint writes per-sample and global tables, then selects samples with viral contigs for annotation. CheckV is skipped for empty viral FASTAs. Provirus coverage comes from the parent Flye contig; its reported length is the extracted region, and it is not marked circular.
