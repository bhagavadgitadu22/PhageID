# VIBRANT and geNomad
rule db_vibrant:
    output: os.path.join(RESULTS_DIR, "logs", "vibrant_db_downloaded.txt")
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs/vibrant_db.log")
    message: "Downloading the VIBRANT database"
    shell: """(date && download-db.sh && echo "Database downloaded" > {output} && date) &> {log}"""

rule vibrant:
    output: os.path.join(RESULTS_DIR, "{sample}", "vibrant", "VIBRANT_filtered_assembly", "VIBRANT_results_filtered_assembly", "VIBRANT_summary_results_filtered_assembly.tsv")
    input: 
        assembly = rules.filtered_assembly_flye.output,
        #db = rules.db_vibrant.output,
        db = "/work/river/Databases/Vibrant/databases"
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_vibrant.log")
    message: "Running VIBRANT"
    shell:
        """(date && VIBRANT_run.py -t {threads} -d {input.db} -i {input.assembly} -folder $(dirname $(dirname $(dirname {output}))) && date) &> {log}"""

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
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_filtered_assembly", "filtered_assembly_summary", "filtered_assembly_virus.fna")
    input: 
        assembly = rules.filtered_assembly_flye.output,
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
        assembly = rules.filtered_assembly_flye.output
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_checkv.log")
    message: "Running the first CheckV per assembly"
    shell:
        """(date && checkv end_to_end -t {threads} -d $(dirname $(dirname {input.db})) {input.assembly} $(dirname {output.checkv_quality}) && date) &> {log}"""

