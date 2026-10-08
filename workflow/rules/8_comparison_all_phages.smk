# Comparison of phages' annotations with lovis4u
rule install_lovis4u_linux:
    output: os.path.join(RESULTS_DIR, "logs", "lovis4u_linux.touch")
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    shell:
        """lovis4u --linux && lovis4u --get-hmms && touch {output}"""

rule lovis4u_pharokka:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "lovis4u_pharokka", "lovis4u.pdf")
    input: 
        linux = rules.install_lovis4u_linux.output,
        gff = active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"))
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "lovis4u_pharokka.log")
    shell:
        r"""(date && 
        mkdir -p $(dirname {output})/gff_input && 
        for gff in {input.gff}; do 
            sample=$(basename $(dirname $(dirname $gff))); 
            ln -sf $(realpath $gff) $(dirname {output})/gff_input/${{sample}}_pharokka.gff; 
        done && 
        lovis4u -gff $(dirname {output})/gff_input --reorient_loci --use-filename-as-id --homology-links --run-hmmscan -o $(dirname {output}) && 
        date) &> {log}"""

rule lovis4u_genotate:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "lovis4u_genotate", "lovis4u.pdf")
    input:
        linux = rules.install_lovis4u_linux.output,
        gff = active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "genotate_pharokka.gff"))
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "lovis4u_genotate.log")
    shell:
        r"""(date && 
        mkdir -p $(dirname {output})/gff_input && 
        for gff in {input.gff}; do 
            sample=$(basename $(dirname $(dirname $(dirname $gff)))); 
            ln -sf $(realpath $gff) $(dirname {output})/gff_input/${{sample}}_genotate_pharokka.gff; 
        done && 
        lovis4u -gff $(dirname {output})/gff_input --reorient_loci --use-filename-as-id --homology-links --run-hmmscan -o $(dirname {output}) && 
        date) &> {log}"""

# theBIGbam for all phages individually
rule thebigbam_mapping:
    output: 
        bam = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam"),
        bai = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam.bai")
    input:
        assembly = rules.fix_circular_viral_contigs_per_sample.output.corrected,
        read1 = rules.preprocess_reads_porechop.output
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_mapping.log")
    message: "Running minimap2 to calculate coverage"
    shell:
        """(date &&
        thebigbam mapping-per-sample -t {threads} -r1 {input.read1:q} -a {input.assembly:q} --circular -o {output.bam:q} &&
        date) &> {log}"""

rule thebigbam_mapping_without_bacteria:
    output: 
        bam = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam"),
        bai = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam.bai")
    input:
        assembly = rules.fix_circular_viral_contigs_per_sample.output.corrected,
        read1 = rules.remove_bacterial_contamination.output
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_mapping_without_bacteria.log")
    message: "Running minimap2 to calculate coverage"
    shell:
        """(date &&
        thebigbam mapping-per-sample -t {threads} -r1 {input.read1:q} -a {input.assembly:q} --circular -o {output.bam:q} &&
        date) &> {log}"""

rule thebigbam_annotations:
    output: os.path.join(RESULTS_DIR, "thebigbam", "viral_annotations_enriched.gff")
    input:
        annotations=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "annotations_on_viruses_enriched", "viral_annotations_enriched.gff")),
        fastas=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"))
    params:
        script=workflow.basedir + "/scripts/filter_dereplicated_annotations.py"
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "combined_viral_annotations.log")
    message: "Combining enriched viral annotations from all samples"
    shell:
        """(date &&
        python {params.script:q} --annotations {input.annotations:q} --fasta {input.fastas:q} --output {output:q} &&
        date) &> {log:q}"""

rule thebigbam_calculate:
    output: os.path.join(RESULTS_DIR, "thebigbam", "thebigbam_ALP.db")
    input:
        annotation = rules.thebigbam_annotations.output,
        bam = active_sample_files(os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam")),
        bam_without_bacteria = active_sample_files(os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam"))
    params:
        bam_dir = os.path.join(RESULTS_DIR, "minimap2", "thebigbam"),
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "thebigbam_calculate.log")
    shell:
        """(date &&
        thebigbam calculate -t {threads} -b {params.bam_dir:q} -g {input.annotation:q} -o {output} --time --blast &&
        date) &> {log}"""
