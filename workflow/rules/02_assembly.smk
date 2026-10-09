# assembling with flye
# using --meta as it is recommended for extrachromosomal elements like phages or plasmids
rule assembly_reads_flye:
    output:
        assembly = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_graph.gfa"),
        info = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_info.txt")
    input: rules.preprocess_reads_porechop.output
    conda: os.path.join(ENV_DIR, "flye.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads_flye.log")
    threads: 8
    shell:
        """(date && flye -t {threads} --meta --nano-raw {input} -o $(dirname {output.assembly}) && date) &> {log}"""

# Single-end short-read assembly, selected automatically when long reads are absent.
rule assembly_reads_spades:
    output: os.path.join(RESULTS_DIR, "{sample}", "spades", "assembly.fasta")
    input: rules.prepare_short_reads.output
    conda: os.path.join(ENV_DIR, "spades.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_spades.log")
    shell:
        """(spades.py -s {input:q} --isolate -t {threads} -o $(dirname {output.assembly:q}) && cp $(dirname {output.assembly:q})/contigs.fasta {output.assembly:q}) > {log:q} 2>&1"""

# Long-read rescue assembly; can also be requested explicitly.
rule autocycler_genome_size:
    input: rules.preprocess_reads_porechop.output
    output: os.path.join(RESULTS_DIR, "{sample}", "autocycler", "genome_size.txt")
    params: genome_size=str(config.get("autocycler_genome_size", "auto"))
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_autocycler_genome_size.log")
    shell:
        """
        if [ {params.genome_size:q} = auto ]; then
            autocycler helper genome_size --reads {input:q} --threads {threads} > {output:q} 2> {log:q}
        else
            printf '%s\n' {params.genome_size:q} > {output:q}
        fi
        """

rule autocycler_subsample:
    input:
        reads=rules.preprocess_reads_porechop.output,
        genome_size=rules.autocycler_genome_size.output
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "autocycler", "read_subsets"))
    params: subset_count=config.get("autocycler_subsets", 4)
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_autocycler_subsample.log")
    shell:
        """autocycler subsample --reads {input.reads:q} --out_dir {output:q} --count {params.subset_count} --genome_size "$(cat {input.genome_size:q})" > {log:q} 2>&1"""

rule autocycler_candidate_assemblies:
    input:
        subsets=rules.autocycler_subsample.output,
        genome_size=rules.autocycler_genome_size.output
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "autocycler", "candidates"))
    params: read_type=config.get("autocycler_read_type", "ont_r10")
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_autocycler_candidates.log")
    shell:
        """
        (
            mkdir -p {output:q}/assemblies {output:q}/logs
            shopt -s nullglob
            subsets=({input.subsets:q}/sample_*.fastq)
            (( ${{#subsets[@]}} > 0 ))
            genome_size=$(cat {input.genome_size:q})
            for assembler in canu flye metamdbg miniasm necat nextdenovo plassembler raven; do
                for reads in "${{subsets[@]}}"; do
                    name="$assembler"_$(basename "$reads" .fastq)
                    prefix={output:q}/"$name"
                    if autocycler helper "$assembler" --reads "$reads" --out_prefix "$prefix" --threads {threads} --genome_size "$genome_size" --read_type {params.read_type:q} > {output:q}/logs/"$name".log 2>&1 && [ -s "$prefix.fasta" ]; then
                        cp "$prefix.fasta" {output:q}/assemblies/"$name".fasta
                    else
                        echo "Assembly failed: $name (see candidate log)."
                    fi
                done
            done
            assemblies=({output:q}/assemblies/*.fasta)
            if (( ${{#assemblies[@]}} < 2 )); then
                echo "Autocycler requires multiple successful input assemblies." >&2
                exit 1
            fi
        ) > {log:q} 2>&1
        """

rule assembly_reads_autocycler:
    input:
        reads=rules.preprocess_reads_porechop.output,
        candidates=rules.autocycler_candidate_assemblies.output
    output: os.path.join(RESULTS_DIR, "{sample}", "autocycler", "assembly.fasta")
    params: work=lambda wc: os.path.join(RESULTS_DIR, wc.sample, "autocycler", "consensus_work")
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_autocycler_consensus.log")
    shell:
        """
        (
            mkdir -p {params.work:q}
            attempt=$(mktemp -d {params.work:q}/attempt_XXXXXX)
            echo "Autocycler working directory: $attempt"
            graph="$attempt/autocycler_out"
            autocycler compress -i {input.candidates:q}/assemblies -a "$graph" -t {threads}
            autocycler cluster -a "$graph"
            shopt -s nullglob
            clusters=("$graph"/clustering/qc_pass/cluster_*)
            if (( ${{#clusters[@]}} == 0 )); then
                echo "No clusters passed Autocycler QC." >&2
                exit 1
            fi
            for cluster in "${{clusters[@]}}"; do
                autocycler trim -c "$cluster" -t {threads}
                autocycler resolve -c "$cluster"
            done
            autocycler combine -a "$graph" -i "$graph"/clustering/qc_pass/cluster_*/5_final.gfa -r {input.reads:q} -t {threads}
            test -s "$graph/consensus_assembly.fasta"
            cp "$graph/consensus_assembly.fasta" {output:q}
        ) > {log:q} 2>&1
        """
