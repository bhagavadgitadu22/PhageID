rule db_pharokka:
    output: directory(os.path.join(RESULTS_DIR, "dbs", "pharokka_v1.4.0_databases"))
    log: os.path.join(RESULTS_DIR, "logs", "pharokka_db.log")
    shell:
        """(date && 
        cd $(dirname {output}) &&
        wget "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz" &&
        tar -xzf pharokka_v1.4.0_databases.tar.gz && 
        date) &> {log}"""

# Protein annotation of a phage
rule pharokka_phage:
    output: os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk")
    input: 
        virus = rules.filtered_assembly.output,
        db = rules.db_pharokka.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka.py --force -t {threads} -d {input.db} -i {input.virus} --dnaapler -o $(dirname {output}) && date) &> {log}"""

rule pharokka_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "pharokka", "plots", "{sample}_annotated_by_pharokka.png")
    input: 
        virus = rules.filtered_assembly.output,
        pharokka_gbk = rules.pharokka_phage.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && pharokka_plotter.py -i {input.virus} -n $(echo {output} | sed 's/.png//') -o $(dirname {input.pharokka_gbk}) && date) &> {log}"""

# Alternative annotation with genotate
rule genotate_calling:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate.faa")
    input: rules.filtered_assembly.output
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
        """(date && sed 's/#//' {input} | sed 's/+//' > {output} && date) &> {log}"""

rule genotate_pharokka:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "pharokka_proteins_full_merged_output.tsv")
    input: 
        virus = rules.genotate_reformatting.output,
        db = rules.db_pharokka.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_annotation.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka_proteins.py --force -t {threads} -d {input.db} -i {input.virus} -o $(dirname {output}) && date) &> {log}"""

rule genotate_pharokka_gff:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "genotate_pharokka.gff")
    input: rules.genotate_pharokka.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_gff_output.log")
    shell:
        """(date && python scripts/reformat_pharokka_genotate.py {input} {output} && date) &> {log}"""

rule genotate_pharokka_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "genotate_pharokka.png")
    input:
        metadata = rules.phage_contig_info.output,
        gff = rules.genotate_pharokka_gff.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_gff_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && python scripts/plotting_pharokka_genotate.py {wildcards.sample} $(tail -n +2 {input.metadata} | cut -f 3) {input.gff} {output} && date) &> {log}"""
