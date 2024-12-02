rule concat_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq")
    input: os.path.join(READS_DIR, "{sample}")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_concat_reads.log")
    shell:
        """cat {input}/*.fastq > {output}"""

rule preprocess_reads_fastp:
    output: 
        fastq = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.{sample}.fastq"),
        html = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.html"),
        json = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.json")
    input: rules.concat_reads.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_fastp.log")
    threads: 4
    shell:
        """(date && fastp -t {threads} -i {input} -o {output.fastq} --html {output.html} --json {output.json} --length_required 1000 --average_qual 20 --cut_front --cut_front_mean_quality 20 && date) &> {log}"""

rule preprocess_reads_fastplong:
    output: 
        fastq = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastplong.{sample}.fastq"),
        html = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastplong.html"),
        json = os.path.join(RESULTS_DIR, "{sample}", "reads", "fastplong.json")
    input: rules.concat_reads.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_fastplong.log")
    threads: 4
    shell:
        """(date && fastplong -t {threads} -i {input} -o {output.fastq} --html {output.html} --json {output.json} --length_required 1000 --mean_qual 20 --cut_front --cut_front_mean_quality 20 && date) &> {log}"""

rule assembly_reads:
    output: 
        assembly = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "assembly.gfa")
    input: rules.preprocess_reads_fastplong.output.fastq
    conda: os.path.join(ENV_DIR, "assembly.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads.log")
    threads: 8
    shell:
        """(date && unicycler -t {threads} -l {input} -o $(dirname {output.assembly}) && date) &> {log}"""

rule depths_per_contig:
    output: 
        depths = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "depths_per_contig.tsv"),
        list_contigs = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "contigs_to_keep.txt")
    input: rules.assembly_reads.output.graph
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_depths_per_contig.log")
    shell:
        """(date && 
        grep -P 'S\t' {input} | cut -f 2,5 | sed 's/dp:f://' > {output.depths} && 
        awk '{{if ($2 >= 0.1 * max) print $1}}' max=$(cut -f 2 {output.depths} | sort -n | tail -1) {output.depths} > {output.list_contigs} &&
        date) &> {log}"""

rule filtered_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "unicycler", "filtered_assembly.fasta")
    input: 
        assembly = rules.assembly_reads.output.assembly,
        list_contigs = rules.depths_per_contig.output.list_contigs
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_depths_per_contig.log")
    shell:
        """(date && seqtk subseq {input.assembly} {input.list_contigs} > {output} && 
        sed -i "s/>/>{wildcards.sample}_/" {output} && date) &> {log}"""

rule quality_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "quast", "report.html")
    input: rules.filtered_assembly.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_quality_assembly.log")
    shell:
        """(date && quast.py {input} -o $(dirname {output}) && date) &> {log}"""

rule multiqc_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "multiqc_report.html")
    input: 
        fastp = rules.preprocess_reads_fastplong.output,
        quast = rules.quality_assembly.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_multiqc_assembly.log")
    message: "Running MultiQC"
    shell:
        """(date && multiqc $(dirname {input.fastp}) $(dirname {input.quast}) --force --config scripts/multiqc_config.yaml -o $(dirname {output}) && date) &> {log}"""
