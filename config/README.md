Edit `config/config.yaml` before running PhageID. The paths in the committed file
describe the original study environment and need to be replaced for another system.

| Setting | What to provide |
| --- | --- |
| `samples_file` | Path to a tab-separated sample sheet, such as `../data/samples.example.tsv`. Relative paths are resolved from `workflow/`. |
| `sample_data_root` | Base directory for relative read and host-reference paths in the sample sheet. |
| `results_dir` | Directory for output files. An absolute path is recommended. |
| `genomad_db`, `checkv_db`, `pharokka_db`, `phold_db`, `plassembler_db`, `checkamg_db`, `vcontact3_db` | Installed reference database directories used by the selected workflow stages. |

The sample sheet has no header and contains four tab-separated columns: sample ID,
long-read FASTQ file or directory, short-read FASTQ file or directory, and host
reference FASTA. Use `NA` for an absent input. Each sample must have at least one
read input. `data/samples.example.tsv` shows the layout. The study sample sheet
in `data/samples.tsv` refers to data that must be downloaded separately.

After configuring paths, activate the environment from `envs/requirements.yaml`
and run from the repository root:

```bash
snakemake --use-conda --cores 8 pipeline_all_phages
```

`corrected_assembly_all_phages` runs assembly and reports only;
`pipeline_comparison_phages` adds the cross-sample comparison. The same commands
also work from `workflow/`, as shown in the main README.
