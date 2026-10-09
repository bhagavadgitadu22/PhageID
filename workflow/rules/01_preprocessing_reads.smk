# Both read types accept a FASTQ file or a directory of FASTQ files.
rule concat_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "{read_type}.{sample}.fastq")
    input: lambda wc: read_input_files((READ_FILES if wc.read_type == "long" else SHORT_READS)[wc.sample])
    wildcard_constraints: read_type="long|short"
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{read_type}_concat_reads.log")
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
    input: lambda wc: rules.concat_reads.output[0].format(sample=wc.sample, read_type="long")
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_porechop.log")
    shell:
        """(date && porechop -t {threads} -i {input} -o {output} --discard_middle --check_reads 10000 --end_size 200 && date) &> {log:q}"""

rule preprocess_reads_fastp:
    output:
        reads=os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.{sample}.fastq"),
        json=os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.{sample}.json"),
        html=os.path.join(RESULTS_DIR, "{sample}", "reads", "fastp.{sample}.html")
    input: lambda wc: rules.concat_reads.output[0].format(sample=wc.sample, read_type="short")
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_fastp.log")
    shell:
        """(date && fastp --thread {threads} -i {input:q} --qualified_quality_phred 20 --length_required 50 \
            -o {output.reads:q} --json {output.json:q} --html {output.html:q} && date) &> {log:q}"""

# remove all reads mapping to the bacterial host genome to keep only phage reads
# full reads remain available for prophage rescue
rule remove_bacterial_contamination:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "cleaned.{sample}.fastq")
    input: 
        reads = sample_reads,
        ref = lambda wc: [HOSTS_LIST[wc.sample]] if HOSTS_LIST[wc.sample] else [],
    params:
        has_host=lambda wc: int(bool(HOSTS_LIST[wc.sample])),
        preset=lambda wc: "map-ont" if READ_FILES[wc.sample] else "sr"
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_remove_bacterial_contamination.log")
    shell:
        """(date
        if (( {params.has_host} )); then
            minimap2 -t {threads} -ax {params.preset} {input.ref:q} {input.reads:q} |
                samtools view -@ {threads} -b -f 4 - | samtools fastq - > {output:q}
        else
            echo "No host genome supplied; retaining all adapter-trimmed reads."
            cp {input.reads:q} {output:q}
        fi
        date) &> {log:q}"""

checkpoint assembly_read_status:
    input:
        reads=sample_reads,
        filtered=rules.remove_bacterial_contamination.output
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "host_filtering_stats.tsv")
    params:
        script=workflow.basedir + "",
        has_host=lambda wc: int(bool(HOSTS_LIST[wc.sample]))
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_host_filtering_stats.log")
    shell:
        """python ./scripts/assembly_read_stats.py filter \
            --all-reads {input.reads:q} --filtered-reads {input.filtered:q} \
            --has-host {params.has_host} --output {output:q} > {log:q} 2>&1"""
