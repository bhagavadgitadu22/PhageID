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

# Termination
rule phageterm_preparation:
    output: directory(os.path.join(ENV_DIR, "PhageTerm", "phagetermvirome-4.3", "phagetermvirome"))
    log: os.path.join(RESULTS_DIR, "logs", "phageterm_preparation.log")
    conda: os.path.join(ENV_DIR, "phageterm.yaml")
    message: "Preparing PhageTermVirome"
    shell:
        """(date && cd $(dirname $(dirname {output})) &&
        pip install poetry &&
        export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True &&
        wget -nc https://files.pythonhosted.org/packages/ef/89/50321c714580c79d431cd9eb12aa62dc49e6f44afbe4e3efae282c9138ff/phagetermvirome-4.3.tar.gz &&
        tar -xzf phagetermvirome-4.3.tar.gz &&
        cd phagetermvirome-4.3 &&
        poetry install &&
        date) &> {log}"""

rule phageterm:
    output: os.path.join(RESULTS_DIR, "{sample}", "phageterm", "Analysis_PhageTerm_report.pdf")
    input:
        reads = rules.preprocess_reads_fastplong.output.fastq,
        virus = rules.filtered_assembly_flye.output,
        phageterm = rules.phageterm_preparation.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phageterm.log")
    conda: os.path.join(ENV_DIR, "phageterm.yaml")
    threads: 8
    shell:
        """(date && export PYTHONPATH={input.phageterm} &&
        cd $(dirname {output}) &&
        phageterm -c {threads} -f {input.reads} -r {input.virus} --report_title Analysis && date) &> {log}"""
