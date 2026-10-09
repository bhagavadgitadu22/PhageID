# dereplication of the viruses
rule combine_all_viruses:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "all_viruses.fna")
    input: active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"))
    log: os.path.join(RESULTS_DIR, "logs", "combine_all_viruses.log")
    message: "Combining all viruses from all samples"
    shell:
        """
        (date && mkdir -p $(dirname {output}) &&
        cat {input} > {output} && date) &> {log}
        """

rule blast_before_dereplication:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "blastn_all_viruses.tsv")
    input: rules.combine_all_viruses.output
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "blastn_all_viruses.log")
    message: "BLAST all-vs-all of the viruses"
    shell:
        """
        (date && blastn -num_threads {threads} -query {input} -subject {input} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out {output} &&
        date) &> {log}
        """

rule ani_for_dereplication:
    output: 
        ani_results = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "ani_all_viruses.tsv"),
        clustering_results = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "clusters_all_viruses.tsv"),
    input:
        blast_results = rules.blast_before_dereplication.output, 
        fna_viruses = rules.combine_all_viruses.output
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "ani_all_viruses.log")
    message: "ANI for dereplication of the viral dataset"
    shell:
        """
        (date && python scripts/ani_calc.py -i {input.blast_results} -o {output.ani_results} &&
        python scripts/ani_clust.py --fna {input.fna_viruses} --ani {output.ani_results} --out {output.clustering_results} --min_ani 95 --min_tcov 85 --min_qcov 0 && date) &> {log}
        """

rule viruses_dereplicated:
    output:
        list_viruses_derep = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "all_viruses_dereplicated.txt"),
        fna_viruses_derep = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "all_viruses_dereplicated.fna"),
    input:
        fna_viruses = rules.combine_all_viruses.output,
        clustering_results = rules.ani_for_dereplication.output.clustering_results,
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "dereplication_all_viruses.log")
    message: "Dereplication of the viral dataset"
    shell:
        """
        (date && cut -f 1 {input.clustering_results} > {output.list_viruses_derep} &&
        seqtk subseq {input.fna_viruses} {output.list_viruses_derep} > {output.fna_viruses_derep} && date) &> {log}
        """

rule dereplication_report:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "dereplication_report.tsv")
    input:
        fasta=rules.combine_all_viruses.output,
        clusters=rules.ani_for_dereplication.output.clustering_results,
        representatives=rules.viruses_dereplicated.output.fna_viruses_derep,
        ani=rules.ani_for_dereplication.output.ani_results
    log: os.path.join(RESULTS_DIR, "logs", "dereplication_report.log")
    shell:
        """python ./scripts/dereplication_report.py --fasta {input.fasta:q} --clusters {input.clusters:q} --representatives {input.representatives:q} --ani {input.ani:q} --output {output:q} > {log:q} 2>&1"""

# Resolve data-dependent sample/representative pairs after clustering.
checkpoint prepare_post_dereplication_mapping:
    output:
        prepared=directory(os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "mapping_inputs"))
    input:
        representatives=rules.viruses_dereplicated.output.fna_viruses_derep,
        clusters=rules.ani_for_dereplication.output.clustering_results,
        sample_fastas=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"))
    params:
        script=workflow.basedir + "/scripts/dereplication_mapping.py",
        sample_args=lambda wildcards, input: [value for sample, fasta in zip(active_viral_samples(wildcards), input.sample_fastas)                                     for value in ("--sample-fasta", sample, str(fasta))]
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "prepare_post_dereplication_mapping.log")
    shell:
        """python {params.script:q} --representatives-fasta {input.representatives:q} --clusters {input.clusters:q} {params.sample_args:q} --output-dir {output.prepared:q} > {log:q} 2>&1"""

def post_dereplication_pairs(wildcards):
    import csv
    prepared = checkpoints.prepare_post_dereplication_mapping.get().output.prepared
    with open(os.path.join(prepared, "sample_representatives.tsv"), newline="") as handle:
        return [(row["sample"], row["representative"]) for row in csv.DictReader(handle, delimiter="\t")]

def post_dereplication_bams(wildcards):
    return [os.path.join(RESULTS_DIR, "minimap2", "thebigbam_post_dereplication", f"{sample}_on_{representative}.bam")
            for sample, representative in post_dereplication_pairs(wildcards)]

def post_dereplication_reference(wildcards):
    if (wildcards.sample, wildcards.representative) not in post_dereplication_pairs(wildcards):
        raise ValueError(f"Sample {wildcards.sample} has no contig in representative cluster {wildcards.representative}")
    prepared = checkpoints.prepare_post_dereplication_mapping.get().output.prepared
    return os.path.join(prepared, "references", wildcards.representative + ".fna")

rule thebigbam_mapping_post_dereplication:
    output:
        bam=os.path.join(RESULTS_DIR, "minimap2", "thebigbam_post_dereplication", "{sample}_on_{representative}.bam")
    input:
        assembly=post_dereplication_reference,
        read1=sample_reads
    wildcard_constraints:
        sample="|".join(__import__("re").escape(sample) for sample in PHAGES_LIST) or "(?!)",
        representative="[^/]+"
    params: mapper=sample_mapper
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_on_{representative}_mapping_post_dereplication.log")
    message: "Mapping {wildcards.sample} reads to representative {wildcards.representative}"
    shell:
        """(date &&
        thebigbam mapping-per-sample -t {threads} -r1 {input.read1:q} -a {input.assembly:q} --mapper {params.mapper:q} --circular -o {output.bam:q} &&
        date) &> {log:q}"""

rule thebigbam_annotations_post_dereplication:
    output: os.path.join(RESULTS_DIR, "thebigbam", "viral_annotations_enriched_post_dereplication.gff")
    input:
        annotations=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "annotations_on_viruses_enriched", "viral_annotations_enriched.gff")),
        representatives=rules.viruses_dereplicated.output.fna_viruses_derep
    params:
        script=workflow.basedir + "/scripts/filter_dereplicated_annotations.py"
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "thebigbam_annotations_post_dereplication.log")
    message: "Combining enriched annotations for dereplicated representatives only"
    shell:
        """(date &&
        python {params.script:q} --annotations {input.annotations:q} --representatives-fasta {input.representatives:q} --output {output:q} &&
        date) &> {log:q}"""

rule thebigbam_calculate_post_dereplication:
    output: os.path.join(RESULTS_DIR, "thebigbam", "thebigbam_ALP_post_dereplication.db")
    input:
        annotation = rules.thebigbam_annotations_post_dereplication.output,
        bam = post_dereplication_bams
    params:
        bam_dir = os.path.join(RESULTS_DIR, "minimap2", "thebigbam_post_dereplication"),
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "thebigbam_calculate_post_dereplication.log")
    shell:
        """(date &&
        thebigbam calculate -t {threads} -b {params.bam_dir:q} -g {input.annotation:q} -o {output} --time --blast &&
        date) &> {log}"""
