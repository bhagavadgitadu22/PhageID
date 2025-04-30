rule db_pharokka:
    output: directory(os.path.join(RESULTS_DIR, "dbs", "pharokka_v1.4.0_databases"))
    log: os.path.join(RESULTS_DIR, "logs", "pharokka_db.log")
    shell:
        """(date && 
        cd $(dirname {output}) &&
        wget "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz" &&
        tar -xzf pharokka_v1.4.0_databases.tar.gz && 
        date) &> {log}"""

rule db_phold:
    output: directory(os.path.join(RESULTS_DIR, "dbs", "phold_databases"))
    log: os.path.join(RESULTS_DIR, "logs", "phold_db.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && 
        cd $(dirname {output}) &&
        phold install -d {output} && 
        date) &> {log}"""

# Protein annotation of a phage
rule pharokka_phage:
    output: 
	    gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
	    gff = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"),
	    dnaapler = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "dnaapler", "dnaapler_reoriented.fasta"),
	    faa = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "phanotate.faa")
    input: 
        virus = rules.filtered_assembly_flye.output,
        db = rules.db_pharokka.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka.py --force -t {threads} -d {input.db} -i {input.virus} --dnaapler -o $(dirname {output.gbk}) && date) &> {log}"""

rule pharokka_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "pharokka", "plots", "{sample}_annotated_by_pharokka.png")
    input: 
        virus = rules.filtered_assembly_flye.output,
        pharokka_gbk = rules.pharokka_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && pharokka_plotter.py -i {input.virus} -n $(echo {output} | sed 's/.png//') -o $(dirname {input.pharokka_gbk}) && date) &> {log}"""

# Annotation with phold in addition
rule phold_phage:
    output: 
	    gbk = os.path.join(RESULTS_DIR, "{sample}", "phold", "phold.gbk")
    input: 
        gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
        db = rules.db_phold.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phold.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && phold run -t {threads} -i {input.gbk} -d --cpu -o $(dirname {output.gbk}) && date) &> {log}"""

rule phold_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "phold", "plots", "{sample}_annotated_by_phold.png")
    input: 
        phold_gbk = rules.phold_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phold_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && phold plot -i {input.phold_gbk} -p $(echo {output} | sed 's/.png//') -o $(dirname {input.phold_gbk}) && date) &> {log}"""
