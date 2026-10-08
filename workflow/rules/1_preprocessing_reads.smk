# concatenate and unzipping source reads into a single fastq file for each sample
rule concat_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq")
    input: lambda wc: READ_FILES[wc.sample],
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_concat_reads.log")
    shell:
        """for f in {input}/*; do
            if [[ "$f" == *.gz ]]; then
                gzip -dc "$f"
            else
                cat "$f"
            fi
        done > {output}"""

# curating the reads to remove adapters and low-quality reads with porechop
rule preprocess_reads_porechop:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "porechop.{sample}.fastq"),
    input: rules.concat_reads.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_porechop.log")
    shell:
        """(date && porechop -t {threads} -i {input} -o {output} --discard_middle --check_reads 10000 --end_size 200 && date) &> {log}"""

# remove all reads mapping to the bacterial host genome to keep only phage reads
# I cannot use in assembly or I would lose prophages though
rule remove_bacterial_contamination:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "cleaned.{sample}.fastq")
    input: 
        reads = rules.preprocess_reads_porechop.output,
        ref = lambda wc: HOSTS_LIST[wc.sample],
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_remove_bacterial_contamination.log")
    shell:
        """(date && minimap2 -t {threads} -ax map-ont {input.ref} {input.reads} | samtools view -@ {threads} -b -f 4 - | samtools fastq - > {output} && date) &> {log}"""
