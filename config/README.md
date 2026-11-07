Edit `config/config.yaml` before running PhageID. The paths in the committed file describe the original study environment and need to be replaced for another system.

| Setting | What to provide |
| --- | --- |
| `samples_file` | Path to a tab-separated sample sheet, such as `../data/samples.example.tsv`. Relative paths are resolved from `workflow/`. |
| `sample_data_root` | Base directory for relative read and host-reference paths in the sample sheet. |
| `results_dir` | Directory for output files. An absolute path is recommended. |
| `genomad_db`, `checkv_db`, `pharokka_db`, `phold_db`, `plassembler_db`, `checkamg_db`, `vcontact3_db` | Installed reference database directories used by the selected workflow stages. |

The sample sheet has no header and contains four tab-separated columns:

- sample ID
  
- long-read FASTQ file or directory
  
- short-read FASTQ file or directory
  
- host reference FASTA.
  

Use `NA` for an absent input. Each sample must have at least one read input. The host genome is optional but can be provided to remove bacterial contamination from the sequencing reads.

`data/samples.example.tsv` shows the layout. The study sample sheet in `data/samples.tsv` refers to the 67 study samples using paths relative to the original ALP collection. Study reads and host reference FASTAs are available under [NCBI BioProject PRJNA1423127](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1423127). You can download and arrange them to match the sample sheet if you want to reproduce the study run.

After configuring paths, activate the environment from `envs/requirements.yaml` and run from the repository root:

```bash
snakemake --use-conda --cores 8 pipeline_all_phages
```

`corrected_assembly_all_phages` runs assembly and reports only; `pipeline_comparison_phages` adds the cross-sample comparison.

If using slurm, run instead:

```bash
SMK_CLUSTER="sbatch --parsable -q {cluster.qos} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time}"
snakemake --use-conda --jobs 200 --local-cores 1 \
  --cluster-config ../config/slurm.yaml \
  --cluster "$SMK_CLUSTER" \
  --cluster-status "python scripts/slurm_cluster_status.py" \
  --keep-going --rerun-incomplete \
  pipeline_all_phages
```

The cluster-status option enables the detection of failed Slurm jobs, including OOM and timeout failures that could otherwise be mistaken for ordinary assembler failures.

More explanations regarding the pipeline can be found in the main [README](https://github.com/bhagavadgitadu22/PhageID).
