# assembling with flye
# using --meta as it is recommended for extrachromosomal elements like phages or plasmids
rule assembly_reads_flye:
    output:
        assembly = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_graph.gfa"),
        info = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_info.txt")
    input: rules.preprocess_reads_porechop.output
    conda: os.path.join(ENV_DIR, "assembly_flye.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads_flye.log")
    threads: 20
    shell:
        """(date && flye -t {threads} --meta --nano-raw {input} -o $(dirname {output.assembly}) && date) &> {log}"""

# running geNomad to check for viral sequences and their taxonomy
rule genomad:
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_filtered_assembly", "filtered_assembly_summary", "filtered_assembly_virus.fna")
    input: 
        assembly = rules.assembly_reads_flye.output.assembly,
        db = config["genomad"]["database"],
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genomad.log")
    message: "Running geNomad"
    shell:
        "(date && genomad end-to-end --threads {threads} --enable-score-calibration --composition virome --max-fdr 0.05 {input.assembly} $(dirname $(dirname {output})) $(dirname {input.db}) && date) &> {log}"

# keeping viral contigs longer than 10 kbp
rule keep_long_viral_contigs:
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_filtered_assembly", "filtered_assembly_summary", "viral_above_10_kbp.fna")
    input: rules.genomad.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_keep_long_viral_contigs.log")
    shell:
        """(date && seqtk seq -L 10000 {input} > {output} && date) &> {log}"""

# renaming contigs with sample name to avoid duplicates in downstream analyses
rule rename_contigs:
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_filtered_assembly", "filtered_assembly_summary", "viral_above_10_kbp_renamed.fna")
    input: rules.keep_long_viral_contigs.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_rename_contigs.log")
    shell:
        """(date && awk -v sample={wildcards.sample} '/^>/ {print ">" sample "_" substr($0, 2); next} {print}' {input} > {output} && date) &> {log}"""

# assessing characteristics of final viral contigs
rule phage_contig_info:
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_stats.tsv")
    input: rules.keep_long_viral_contigs.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_stats.log")
    shell:
        """(date && ./scripts/contig_info.sh -m 1000 -t {input} > {output} && date) &> {log}"""
