# geNomad
rule db_genomad:
    output: os.path.join("/work/river/Databases", "genomad_db", "genomad_marker_metadata.tsv")
    log: os.path.join(RESULTS_DIR, "logs", "genomad_db.log")
    conda: os.path.join(ENV_DIR, "viral_taxonomy.yaml")
    message: "Downloading the geNomad database"
    shell:
        """
        (date && cd $(dirname $(dirname {output})) &&
        wget -nc https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz && 
        tar --skip-old-files -zxvf genomad_db_v1.9.tar.gz && date) &> {log}
        """

rule genomad:
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_filtered_assembly", "autoblast_corrected_breaking_terminal_repeats_summary", "autoblast_corrected_breaking_terminal_repeats_virus.fna")
    input: 
        assembly = rules.break_terminal_repeats_with_autoblast.output.corrected,
        db = rules.db_genomad.output,
    conda: os.path.join(ENV_DIR, "viral_taxonomy.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genomad.log")
    message: "Running geNomad"
    shell:
        "(date && genomad end-to-end --threads {threads} --restart --enable-score-calibration --composition metagenome --max-fdr 0.05 {input.assembly} $(dirname $(dirname {output})) $(dirname {input.db}) && date) &> {log}"

# CheckV to assess completeness
rule db_checkv:
    output: os.path.join("/work/river/Databases", "checkv-db-v1.5", "genome_db", "checkv_reps.dmnd")
    log: os.path.join(RESULTS_DIR, "logs", "checkv_db.log")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
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
        assembly = rules.break_terminal_repeats_with_autoblast.output.corrected
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_checkv.log")
    message: "Running the first CheckV per assembly"
    shell:
        """(date && checkv end_to_end -t {threads} -d $(dirname $(dirname {input.db})) {input.assembly} $(dirname {output.checkv_quality}) && date) &> {log}"""
