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
        """(date && pharokka_plotter.py -i {input.virus} -n $(echo {output} | sed 's/.png//') --truncate 30 -o $(dirname {input.pharokka_gbk}) && date) &> {log}"""
