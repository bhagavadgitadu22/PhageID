#rule db_pharokka:
#     output: directory(os.path.join(RESULTS_DIR, "dbs", "pharokka_v1.4.0_databases"))
#    log: os.path.join(RESULTS_DIR, "logs", "pharokka_db.log")
#    shell:
#        """(date && 
#        cd $(dirname {output}) &&
#        wget "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz" &&
#        tar -xzf pharokka_v1.4.0_databases.tar.gz && 
#        date) &> {log}"""

rule db_phold:
    output: directory(os.path.join("/work/river/Databases", "phold_db"))
    log: os.path.join(RESULTS_DIR, "logs", "phold_db.log")
    conda: os.path.join(ENV_DIR, "phold.yaml")
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
	    faa = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "phanotate.faa")
    input: 
        virus = rules.correct_phage_with_autoblast.output.corrected,
        #db = rules.db_pharokka.output
        db = "/work/river/Databases/pharokka_db"
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka.py --force -t {threads} -d {input.db} -i {input.virus} -o $(dirname {output.gbk}) && date) &> {log}"""

rule pharokka_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "pharokka", "plots", "{sample}_annotated_by_pharokka.png")
    input: 
        virus = rules.correct_phage_with_autoblast.output.corrected,
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
    conda: os.path.join(ENV_DIR, "phold.yaml")
    threads: 8
    shell:
        """(date && phold run -t {threads} --force -i {input.gbk} -d {input.db} --cpu -o $(dirname {output.gbk}) && date) &> {log}"""

rule phold_plot:
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "phold", "plots"))
    input: 
        phold_gbk = rules.phold_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phold_plot.log")
    conda: os.path.join(ENV_DIR, "phold.yaml")
    shell:
        """(date && phold plot --force -i {input.phold_gbk} -o {output} && date) &> {log}"""
