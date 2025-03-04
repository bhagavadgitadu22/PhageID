rule concat_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq")
    input: os.path.join(READS_DIR, "{sample}")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_concat_reads.log")
    shell:
        """cat {input}/*.fastq > {output}"""

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
        """(date && fastplong -t {threads} -i {input} -o {output.fastq} --html {output.html} --json {output.json} --length_required 1000 && date) &> {log}"""

# Assembling
rule assembly_reads_unicycler:
    output: 
        assembly = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "assembly.gfa")
    input: rules.preprocess_reads_fastplong.output.fastq
    conda: os.path.join(ENV_DIR, "assembly.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads_unicycler.log")
    threads: 8
    shell:
        """(date && unicycler -t {threads} -l {input} -o $(dirname {output.assembly}) && date) &> {log}"""

rule assembly_reads_flye:
    output:
        assembly = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_graph.gfa"),
        info = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_info.txt")
    input: rules.preprocess_reads_fastplong.output.fastq
    conda: os.path.join(ENV_DIR, "assembly_flye.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads_flye.log")
    threads: 20
    shell:
        """(date && flye -t {threads} --meta --nano-raw {input} -o $(dirname {output.assembly}) && date) &> {log}"""

# Filtering best contigs
rule depths_per_contig_unicycler:
    output: 
        depths = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "depths_per_contig.tsv"),
        list_contigs = os.path.join(RESULTS_DIR, "{sample}", "unicycler", "contigs_to_keep.txt")
    input: rules.assembly_reads_unicycler.output.graph
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_depths_per_contig_unicycler.log")
    shell:
        """(date && 
        grep -P 'S\t' {input} | cut -f 2,5 | sed 's/dp:f://' > {output.depths} && 
        awk '{{if ($2 >= 0.1 * max) print $1}}' max=$(cut -f 2 {output.depths} | sort -n | tail -1) {output.depths} > {output.list_contigs} &&
        date) &> {log}"""

rule filtered_assembly_unicycler:
    output: os.path.join(RESULTS_DIR, "{sample}", "unicycler", "filtered_assembly.fasta")
    input: 
        assembly = rules.assembly_reads_unicycler.output.assembly,
        list_contigs = rules.depths_per_contig_unicycler.output.list_contigs
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_filtered_assembly_unicycler.log")
    shell:
        """(date && seqtk subseq {input.assembly} {input.list_contigs} > {output} && 
        sed -i "s/>/>{wildcards.sample}_/" {output} && date) &> {log}"""

rule depths_per_contig_flye:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "contigs_to_keep.txt")
    input: rules.assembly_reads_flye.output.info
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_depths_per_contig_flye.log")
    shell:
        """(date &&
        # Count the number of contigs
        contig_count=$(awk 'NR>1 {{print $1}}' {input} | wc -l)

        # Extract the sequence names of contigs that circularize
        circular_contigs=$(awk 'NR>1 && $4 == "Y" {{print $1}}' {input})
        circular_count=$(echo "$circular_contigs" | wc -l)

        # Check if there is only one contig or only one circular contig
        if [ "$contig_count" -eq 1 ]; then
            # Extract the sequence name of the single contig
            awk 'NR>1 {{print $1}}' {input} > {output}
        elif [ "$circular_count" -eq 1 ]; then
            # Write the number of the single circular contig
            echo "$circular_contigs" > {output} 
        fi &&
        date) &> {log}"""

rule filtered_assembly_flye:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "filtered_assembly.fasta")
    input:
        assembly = rules.assembly_reads_flye.output.assembly,
        list_contigs = rules.depths_per_contig_flye.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_filtered_assembly_flye.log")
    shell:
        """(date && seqtk subseq {input.assembly} {input.list_contigs} > {output} &&
        sed -i "s/>/>{wildcards.sample}_/" {output} && date) &> {log}"""

# Assessing final contigs
rule quality_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "quast", "report.html")
    input: rules.filtered_assembly_flye.output
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

rule phage_contig_info:
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_stats.tsv")
    input: rules.filtered_assembly_flye.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_stats.log")
    shell:
        """(date && ./scripts/contig_info.sh -m 1000 -t {input} > {output} && date) &> {log}"""
