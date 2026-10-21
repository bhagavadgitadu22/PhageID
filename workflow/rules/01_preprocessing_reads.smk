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

# trying to see if a peak of sizes indicates expected number of viruses and their sizes
# that only works for long reads though
rule extract_reads_longer_than_20kbp:
    output: os.path.join(RESULTS_DIR, "{sample}", "long_reads", "longer_than_20kbp.{sample}.fasta")
    input: rules.preprocess_reads_porechop.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_extract_reads_longer_than_20kbp.log")
    shell:
        """(date && seqtk seq -a -L 20000 {input:q} > {output:q} && date) &> {log:q}"""

rule genomad_long_reads:
    output: os.path.join(RESULTS_DIR, "{sample}", "long_reads", "genomad_longer_than_20kbp.{sample}", "longer_than_20kbp.{sample}_summary", "longer_than_20kbp.{sample}_virus_summary.tsv")
    input: 
        reads = rules.extract_reads_longer_than_20kbp.output,
        db = "/work/river/Databases/genomad_db/genomad_marker_metadata.tsv"
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    threads: 10
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genomad_long_reads.log")
    shell:
        """(
            date
            if grep -q '^>' {input.reads:q}; then
                genomad end-to-end --threads {threads} --enable-score-calibration --composition virome --max-fdr 0.05 {input.reads:q} $(dirname $(dirname {output:q})) $(dirname {input.db:q})
            else
                echo "No reads at least 20 kbp; skipping geNomad."
                mkdir -p $(dirname {output:q})
                printf 'seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n' > {output:q}
            fi
            date
        ) &> {log:q}"""

rule histogram_long_viral_reads:
    output: os.path.join(RESULTS_DIR, "reports", "viral_long_reads_histograms", "{sample}_viral_long_reads_histogram.txt")
    input: rules.genomad_long_reads.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_histogram_long_viral_reads.log")
    shell:
        """(date &&
        awk -F'\t' 'NR>1 {{v=int($2/1000)*1000; count[v]++; if(v>max) max=v}}
        END {{
            for(i=0; i<=max; i+=1000) {{
                printf "%6d-%6d | ", i, i+999
                for(j=0; j<count[i]; j++) printf "*"
                printf " (%d)\\n", count[i]
            }}
        }}' {input:q} > {output:q} &&
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

rule all_reads_stats:
    input: 
        expand(os.path.join(RESULTS_DIR, "{sample}", "reads", "host_filtering_stats.tsv"), sample=PHAGES_LIST),
        expand(os.path.join(RESULTS_DIR, "reports", "viral_long_reads_histograms", "{sample}_viral_long_reads_histogram.txt"), sample=[sample for sample in PHAGES_LIST if READ_FILES[sample]])
