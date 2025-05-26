# Alternative annotation with genotate
rule genotate_calling:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate.faa")
    input: rules.pharokka_phage.output.dnaapler
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

# Annotation of the genotate CDS
rule genotate_phold_predict:
    output:
        prediction_dir = directory(os.path.join(RESULTS_DIR, "{sample}", "genotate", "prediction_phold"))
    input:
        genes = rules.genotate_reformatting.output,
        db = rules.db_phold.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_annotation.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && phold predict --force -t {threads} -d {input.db} -i {input.genes} -o {output.prediction_dir} && date) &> {log}"""

rule genotate_phold_compare:
    output: 
        annotation_dir = directory(os.path.join(RESULTS_DIR, "{sample}", "genotate", "annotation_phold"))
    input: 
        genes = rules.genotate_reformatting.output,
        prediction_dir = os.path.join(RESULTS_DIR, "{sample}", "genotate", "prediction_phold"),
        db = rules.db_phold.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_annotation.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && phold proteins-compare --force -t {threads} -d {input.db} -i {input.genes} --predictions_dir {output.prediction_dir} -o {output.annotation_dir} && date) &> {log}"""
