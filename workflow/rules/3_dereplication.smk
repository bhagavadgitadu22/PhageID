# Copy all fasta and gff from phages of a host in a common directory
rule copy_files_genus:
    output: 
        fasta = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "unicycler_sequences")),
        gff = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "pharokka_annotations")),
        pharokka_plot = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "pharokka_plots"))
    input: 
        # Match all directories starting with genus
        lambda wildcards: sorted(
            [os.path.join(RESULTS_DIR, d, "unicycler", "filtered_assembly.fasta") 
            for d in os.listdir(RESULTS_DIR) 
            if d.startswith(wildcards.genus + "_") 
            and os.path.exists(os.path.join(RESULTS_DIR, d, "unicycler", "filtered_assembly.fasta"))]
        )
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_copy_files.log")
    shell:
        """(date && mkdir -p $(dirname {output.fasta} | xargs dirname) && mkdir -p $(dirname {output.fasta}) &&
            mkdir -p {output.fasta} && mkdir -p {output.gff} && mkdir -p {output.pharokka_plot} &&
            for f in {input}; do 
                name=$(dirname $f | xargs dirname | xargs basename);
                ln -s $f {output.fasta}/unicycler_${{name}}.fasta;
                gff_file=$(echo $f | sed 's/unicycler\/filtered_assembly.fasta/pharokka\/pharokka.gff/');
                plot_file=$(echo $f | sed 's/unicycler\/filtered_assembly.fasta/pharokka\/plots\/${{name}}_annotated_by_pharokka.png/');
                ln -s $gff_file {output.gff}/pharokka_${{name}}.gff;
                ln -s $plot_file {output.pharokka_plot}/${{name}}_annotated_by_pharokka.png;
            done && date) &> {log}"""

# Comparison of phages' sequences relative to this host
rule fast_ani:
    output: 
        list = os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "list_sequences.txt"),
        fastani = os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "fastani.out")
    input: rules.copy_files_genus.output.fasta
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_fastani.log")
    shell:
        """(date && find {input} -name "*.fasta" > {output.list} &&
        fastANI --rl {output.list} --ql {output.list} --matrix -o {output.fastani} && date) &> {log}"""

# Comparison of phages' annotations relative to this host
rule install_lovis4u_linux:
    output: os.path.join(RESULTS_DIR, "logs", "lovis4u_linux.touch")
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    shell:
        """lovis4u --linux && lovis4u --get-hmms && touch {output}"""

rule lovis4u:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "lovis4u", "lovis4u.pdf")
    input: 
        linux = rules.install_lovis4u_linux.output,
        gff = rules.copy_files_genus.output.gff
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_lovis4u.log")
    shell:
        """(date && lovis4u -gff {input.gff} --reorient_loci --use-filename-as-id --homology-links --set-category-colour --run-hmmscan -o $(dirname {output}) && date) &> {log}"""
