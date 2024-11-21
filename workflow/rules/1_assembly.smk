rule concat_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq")
    input: os.path.join(READS_DIR, "{sample}")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_concat_reads.log")
    shell:
        """(date && gzip -d {input}/*.fastq.gz && cat {input}/*.fastq > {output} && date) &> {log}"""

rule preprocess_reads:
    output: 
        fastq = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.{sample}.fastq"),
        html = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.html"),
        json = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.json")
    input: rules.concat_reads.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads.log")
    shell:
        """(date && fastp -i {input} -o {output.fastq} --html {output.html} --json {output.json} --length_required 1000 --average_qual 20 --cut_front --cut_front_mean_quality 20 && date) &> {log}"""

rule assembly_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "unicycler", "assembly.fasta")
    input: rules.preprocess_reads.output.fastq
    conda: os.path.join(ENV_DIR, "assembly.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads.log")
    shell:
        """(date && unicycler -l {input} -o $(dirname {output}) && date) &> {log}"""

rule quality_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "unicycler", "quast", "report.html")
    input: rules.assembly_reads.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_quality_assembly.log")
    shell:
        """(date && quast.py {input} -o $(dirname {output}) && date) &> {log}"""

rule multiqc_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "multiqc_report.html")
    input: 
        fastp = rules.preprocess_reads.output,
        quast = rules.quality_assembly.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_multiqc_assembly.log")
    message: "Running MultiQC"
    shell:
        """(date && multiqc $(dirname {input.fastp}) $(dirname {input.quast}) --force --config scripts/multiqc_config.yaml -o $(dirname {output}) && date) &> {log}"""
