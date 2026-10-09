# Both read types accept a FASTQ file or a directory of FASTQ files.
rule concat_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq")
    input: lambda wc: read_input_files(READ_FILES[wc.sample])
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_concat_reads.log")
    shell: 
        """(
            for f in {input:q}; do
                if [[ "$f" == *.gz || "$f" == *.GZ ]]; then
                    gzip -dc "$f"
                else
                    cat "$f"
                fi
            done > {output:q}
            if [ ! -s {output:q} ]; then
                echo "Read input is empty." >&2
                exit 1
            fi
        ) 2> {log:q}"""

# curating the reads to remove adapters and low-quality reads with porechop
rule preprocess_reads_porechop:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "porechop.{sample}.fastq"),
    input: rules.concat_reads.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_porechop.log")
    shell:
        """(date && porechop -t {threads} -i {input} -o {output} --discard_middle --check_reads 10000 --end_size 200 && date) &> {log}"""

# TruSeq input is single-end; combine all files into one read library.
rule prepare_short_reads:
    input: lambda wc: read_input_files(SHORT_READS[wc.sample])
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "short.{sample}.fastq")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_short_reads.log")
    shell:
        """(
            for f in {input:q}; do
                if [[ "$f" == *.gz || "$f" == *.GZ ]]; then
                    gzip -dc "$f"
                else
                    cat "$f"
                fi
            done > {output:q}
            if [ ! -s {output:q} ]; then
                echo "Read input is empty." >&2
                exit 1
            fi
        ) 2> {log:q}"""
