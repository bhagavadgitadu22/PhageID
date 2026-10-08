# CheckV to assess completeness
rule db_checkv:
    output: os.path.join(RESULTS_DIR, "dbs", "checkv-db-v1.5", "genome_db", "checkv_reps.dmnd")
    log: os.path.join(RESULTS_DIR, "logs", "checkv_db.log")
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    message: "Downloading the CheckV database"
    shell:
        """
        (date && mkdir -p $(dirname $(dirname {output})) && cd $(dirname $(dirname $(dirname {output}))) &&
        wget -nc https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz &&
        tar --skip-old-files -zxvf checkv-db-v1.5.tar.gz && 
        cd checkv-db-v1.5/genome_db && diamond makedb --in checkv_reps.faa --db checkv_reps && date) &> {log}
        """

rule checkv:
    output: 
        checkv_quality = os.path.join(RESULTS_DIR, "{sample}", "checkv", "quality_summary.tsv"),
    input: 
        db = rules.db_checkv.output,
        assembly = rules.filtered_assembly_flye.output
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_checkv.log")
    message: "Running the first CheckV per assembly"
    shell:
        """(date && checkv end_to_end -t {threads} -d $(dirname $(dirname {input.db})) {input.assembly} $(dirname {output.checkv_quality}) && date) &> {log}"""
