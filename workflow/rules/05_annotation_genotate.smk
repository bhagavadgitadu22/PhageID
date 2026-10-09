# Alternative annotation with genotate
rule genotate_calling:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate.faa")
    input: rules.select_viral_assembly.output.corrected
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_calling.log")
    conda: os.path.join(ENV_DIR, "genotate.yaml")
    threads: 4
    shell:
        """(date && genotate.py {input} -o {output} --format faa && date) &> {log}"""

rule genotate_reformatting:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_reformatted.faa")
    input: rules.genotate_calling.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_reformatting.log")
    shell:
        """(date && sed 's/#//g' {input} | sed 's/+//g' > {output} && date) &> {log}"""

rule genotate_pharokka:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "pharokka_proteins_full_merged_output.tsv")
    input:
        genes = rules.genotate_reformatting.output,
        db = "/work/river/Databases/pharokka_db"
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_annotation.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka proteins --force -t {threads} -d {input.db} -i {input.genes} -o $(dirname {output}) && date) &> {log}"""

rule genotate_pharokka_gff:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "genotate_pharokka.gff")
    input:
        genes = rules.genotate_pharokka.output,
        fasta = rules.select_viral_assembly.output.corrected
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_gff_output.log")
    shell:
        """(date && python ./scripts/reformat_pharokka_genotate.py --genes {input.genes:q} --fasta {input.fasta:q} --output {output:q} && date) &> {log}"""

rule genotate_pharokka_plot:
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "plots"))
    input:
        gff_genotate = rules.genotate_pharokka_gff.output,
        gff_pharokka = rules.pharokka_phage.output.gff,
        fasta = rules.select_viral_assembly.output.corrected
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_gff_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && python ./scripts/plotting_pharokka_genotate.py --fasta {input.fasta:q} --genotate-gff {input.gff_genotate:q} --pharokka-gff {input.gff_pharokka:q} --output-dir {output:q} && date) &> {log}"""
