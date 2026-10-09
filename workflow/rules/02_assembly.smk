# External sampling preserves --meta and leaves full reads available elsewhere.
rule subsample_reads_flye:
    input: assembly_reads
    output:
        reads=os.path.join(RESULTS_DIR, "{sample}", "reads", "flye.{sample}.{readset}.fastq"),
        report=os.path.join(RESULTS_DIR, "{sample}", "reads", "flye_subsampling.{readset}.tsv")
    params:
        script=workflow.basedir + "/scripts/subsample_flye_reads.py",
        genome_size=100000,
        coverage=1000,
        seed=42
    conda: os.path.join(ENV_DIR, "flye.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_flye_subsampling.log")
    shell:
        """python {params.script:q} --reads {input:q} --genome-size {params.genome_size} --coverage {params.coverage} \
            --seed {params.seed} --output {output.reads:q} --report {output.report:q} > {log:q} 2>&1"""

# assembling with flye
# using --meta as it is recommended for extrachromosomal elements like phages or plasmids
rule assembly_reads_flye:
    output:
        assembly = os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "flye", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "flye", "assembly_graph.gfa"),
        info = os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "flye", "assembly_info.txt"),
        status = os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "flye", "assembly_status.txt")
    input: rules.subsample_reads_flye.output.reads
    conda: os.path.join(ENV_DIR, "flye.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_assembly_reads_flye.log")
    threads: 8
    shell:
        """python ./scripts/run_assembly_attempt.py --assembly {output.assembly:q} --status {output.status:q} --log {log:q} \
            --extra-output {output.graph:q} --extra-output {output.info:q} -- \
            flye -t {threads} --meta --nano-raw {input:q} -o $(dirname {output.assembly:q})"""

# Single-end short-read assembly, selected automatically when long reads are absent.
rule assembly_reads_spades:
    output:
        assembly=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "spades", "assembly.fasta"),
        status=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "spades", "assembly_status.txt")
    input: assembly_reads
    conda: os.path.join(ENV_DIR, "spades.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_spades.log")
    shell:
        """python ./scripts/run_assembly_attempt.py --assembly {output.assembly:q} --status {output.status:q} --log {log:q} -- \
            bash -euo pipefail -c 'spades.py -s "$1" --isolate -t "$2" -o "$3" && cp "$3/contigs.fasta" "$4"' \
            -- {input:q} {threads} $(dirname {output.assembly:q}) {output.assembly:q}"""
