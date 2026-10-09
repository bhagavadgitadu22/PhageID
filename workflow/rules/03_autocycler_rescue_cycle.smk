# Long-read rescue assembly
rule autocycler_genome_size:
    input: assembly_reads
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "autocycler", "genome_size.txt")
    params: genome_size=str(config.get("autocycler_genome_size", "auto"))
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_autocycler_genome_size.log")
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
        reads=assembly_reads,
        genome_size=rules.autocycler_genome_size.output
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "autocycler", "read_subsets"))
    params: subset_count=4
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_autocycler_subsample.log")
    shell:
        """autocycler subsample --reads {input.reads:q} --out_dir {output:q} --count {params.subset_count} --genome_size "$(cat {input.genome_size:q})" > {log:q} 2>&1"""

def autocycler_candidate_jobs(wc):
    subsets = int(rules.autocycler_subsample.params.subset_count)
    if subsets < 1:
        raise ValueError("subset_count must be at least 1")
    return expand(rules.autocycler_candidate_assembly.output[0], sample=wc.sample, readset=wc.readset,
                  candidate_assembler=["canu", "flye", "metamdbg", "miniasm", "necat", "nextdenovo", "plassembler", "raven"],
                  subset=[f"{i:02d}" for i in range(1, subsets + 1)])

# Each job runs one assembler on one subset. Failed attempts remain reviewable.
rule autocycler_candidate_assembly:
    input:
        subsets=rules.autocycler_subsample.output[0],
        genome_size=rules.autocycler_genome_size.output
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "autocycler", "candidate_jobs", "{candidate_assembler}", "sample_{subset}"))
    wildcard_constraints:
        candidate_assembler="canu|flye|metamdbg|miniasm|necat|nextdenovo|plassembler|raven",
        subset="[0-9]+"
    params:
        read_type="ont_r10",
        plassembler_db="/work/river/Databases/plassembler_db",
        reads=lambda wc, input: os.path.join(input.subsets, f"sample_{wc.subset}.fastq")
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_autocycler_{candidate_assembler}_{subset}.log")
    shell:
        """
        export PLASSEMBLER_DB={params.plassembler_db:q}
        test -s {params.reads:q}
        mkdir -p {output:q}
        extra_args=()
        if [ {wildcards.candidate_assembler:q} = plassembler ]; then
            extra_args=(--args --no_chromosome)
        fi
        if autocycler helper {wildcards.candidate_assembler:q} --threads {threads} --reads {params.reads:q} --out_prefix {output:q}/assembly --genome_size "$(cat {input.genome_size:q})" --read_type {params.read_type:q} "${{extra_args[@]}}" > {log:q} 2>&1 && [ -s {output:q}/assembly.fasta ]; then
            printf 'success\n' > {output:q}/status.txt
        else
            rm -f {output:q}/assembly.fasta
            printf 'failed\n' > {output:q}/status.txt
            echo "Candidate failed; see {log}." >&2
        fi
        """

# Wait for all attempts, then collect only successful FASTAs for consensus.
rule autocycler_candidate_assemblies:
    input: autocycler_candidate_jobs
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "autocycler", "candidates"))
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_autocycler_candidates.log")
    shell:
        """python ./scripts/collect_autocycler_candidates.py --candidates {input:q} --output {output:q} > {log:q} 2>&1"""

rule assembly_reads_autocycler:
    input:
        reads=assembly_reads,
        candidates=rules.autocycler_candidate_assemblies.output
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "autocycler", "assembly.fasta")
    params: work=lambda wc: os.path.join(RESULTS_DIR, wc.sample, "assembly_attempts", wc.readset, "autocycler", "consensus_work")
    conda: os.path.join(ENV_DIR, "autocycler.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_autocycler_consensus.log")
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
