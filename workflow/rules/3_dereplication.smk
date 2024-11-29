GENUS_PATTERN = f"{RESULTS_DIR}/{{genus}}_{{name}}"

rule copy_gff_files:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "gffs")
    input: expand(os.path.join(RESULTS_DIR, "{genus_phages}", "pharokka", "pharokka.gff"), genus_phages=glob_wildcards(GENUS_PATTERN).genus)
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_copy_gff_files.log")
    shell:
        """(date && for f in {input}; do cp $f {output}; done && date) &> {log}"""

rule lovis4u:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "lovis4u")
    input: rules.copy_gff_files.output
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_lovis4u.log")
    shell:
        """(date && lovis4u -gff {input} --reorient_loci --homology-links --set-category-colour -o {output} && date) &> {log}"""
