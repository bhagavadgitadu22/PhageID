# Comparison of phages' annotations with lovis4u
rule install_lovis4u_linux:
    output: os.path.join(RESULTS_DIR, "logs", "lovis4u_linux_pyhmmer_0.11.4.touch")
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    shell:
        """lovis4u --linux && lovis4u --get-hmms && touch {output}"""

rule lovis4u_pharokka:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "lovis4u_pharokka", "lovis4u.pdf")
    input: 
        linux = rules.install_lovis4u_linux.output,
        gff = active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff")),
        splitter = "scripts/prepare_lovis4u_gff.py"
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "lovis4u_pharokka.log")
    shell:
        """(date &&
        python {input.splitter:q} --output-dir $(dirname {output:q})/gff_input {input.gff:q} &&
        lovis4u -gff $(dirname {output:q})/gff_input --reorient_loci --use-filename-as-id --add-locus-id-prefix --homology-links --run-hmmscan --clust_loci-off --debug -o $(dirname {output:q}) &&
        date) &> {log:q}"""

# One annotation/assembly file per sample defines one "viral MAG" in the first database.
rule thebigbam_annotations:
    output: directory(os.path.join(RESULTS_DIR, "thebigbam", "mag_inputs"))
    input:
        annotations=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "annotations_on_viruses_enriched", "viral_annotations_enriched.gff")),
        fastas=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"))
    log: os.path.join(RESULTS_DIR, "logs", "combined_viral_annotations.log")
    message: "Preparing one viral MAG per sample"
    shell:
        """
        (
            mkdir -p {output:q}/annotations {output:q}/assemblies
            for gff in {input.annotations:q}; do
                sample=$(basename "$(dirname "$(dirname "$gff")")")
                cp "$gff" {output:q}/annotations/"$sample".gff
            done
            for fasta in {input.fastas:q}; do
                sample=$(basename "$(dirname "$(dirname "$fasta")")")
                cp "$fasta" {output:q}/assemblies/"$sample".fasta
            done
        ) > {log:q} 2>&1
        """

rule thebigbam_calculate:
    output: os.path.join(RESULTS_DIR, "thebigbam", "thebigbam_ALP.db")
    input:
        annotation = rules.thebigbam_annotations.output,
        bam = active_sample_files(os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam")),
        bam_without_bacteria = active_sample_files(os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam")),
        bai = active_sample_files(os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam.bai")),
        bai_without_bacteria = active_sample_files(os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam.bai"))
    params:
        bam_dir = os.path.join(RESULTS_DIR, "minimap2", "thebigbam"),
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "thebigbam_calculate.log")
    shell:
        """(date &&
        thebigbam calculate -t {threads} -b {params.bam_dir:q} \
            -g {input.annotation:q}/annotations -a {input.annotation:q}/assemblies \
            --view mag -m coverage,misalignment,longreads,phagetermini \
            -o {output:q} --diagnose --blast all &&
        date) &> {log}"""
