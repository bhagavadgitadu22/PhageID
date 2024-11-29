# Copy all fasta and gff from phages of a host in a common directory
rule copy_fasta_files_genus:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "sequences")
    input: 
        # Match all directories starting with genus
        lambda wildcards: sorted(
            [os.path.join(RESULTS_DIR, d, "unicycler", "filtered_assembly.fasta") for d in os.listdir(RESULTS_DIR) if d.startswith(wildcards.genus + "_")]
        )
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_copy_fasta_files.log")
    shell:
        """(date && for f in {input}; do 
                name=$(sed 's/unicycler\/filtered_assembly.fasta//' $f | basename);
                ln -s $f {output}/unicycler_${name}.fasta; 
            done && date) &> {log}"""

rule copy_gff_files_genus:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "annotations")
    input: 
        # Match all directories starting with genus
        lambda wildcards: sorted(
            [os.path.join(RESULTS_DIR, d, "pharokka", "pharokka.gff") for d in os.listdir(RESULTS_DIR) if d.startswith(wildcards.genus + "_")]
        )
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_copy_gff_files.log")
    shell:
        """(date && for f in {input}; do 
                name=$(sed 's/pharokka\/pharokka.gff//' $f | basename);
                ln -s $f {output}/pharokka_${name}.gff; 
            done && date) &> {log}"""

# Comparison of phages' sequences relative to this host
rule fast_ani:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}_fastani.out")
    input: rules.copy_fasta_files_genus.output
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_lovis4u.log")
    shell:
        """(date && fastANI -r {input} -q {input} --matrix -o {output} && date) &> {log}"""

# Comparison of phages' annotations relative to this host
rule lovis4u:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "lovis4u")
    input: rules.copy_gff_files_genus.output
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_lovis4u.log")
    shell:
        """(date && lovis4u -gff {input} --reorient_loci --homology-links --set-category-colour -o {output} && date) &> {log}"""
