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
    threads: 20
    shell:
        """(date && flye -t {threads} --meta --nano-raw {input} -o $(dirname {output.assembly}) && date) &> {log}"""

# running geNomad to check for viral sequences and their taxonomy
rule genomad:
    output:
        fasta=os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus.fna"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus_summary.tsv")
    input: 
        assembly = rules.assembly_reads_flye.output.assembly,
        db = "/work/river/Databases/genomad_db/genomad_marker_metadata.tsv"
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genomad.log")
    message: "Running geNomad"
    shell:
        "(date && genomad end-to-end --threads {threads} --enable-score-calibration --composition virome --max-fdr 0.05 {input.assembly} $(dirname $(dirname {output.fasta})) $(dirname {input.db}) && date) &> {log}"

# keeping viral contigs longer than 2 kbp
rule keep_long_viral_contigs:
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "viral_above_2_kbp.fna")
    input: rules.genomad.output.fasta
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_keep_long_viral_contigs.log")
    shell:
        """(date && seqtk seq -L 2000 {input} > {output} && date) &> {log}"""

# renaming contigs with sample name to avoid duplicates in downstream analyses
rule rename_contigs:
    output: os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "viral_above_2_kbp_renamed.fna")
    input: rules.keep_long_viral_contigs.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_rename_contigs.log")
    shell:
        """(date && awk -v sample={wildcards.sample} '/^>/ {{print ">" sample "_" substr($0, 2); next}} {{print}}' {input} > {output} && date) &> {log}"""

# autoblast to detect potential duplications of the phage (concatemers)
rule fix_circular_viral_contigs_per_sample:
    output:
        corrected=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"),
        concatemer_report=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "breaking_concatemers_report.csv"),
        dtr_report=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "breaking_dtr_report.csv")
    input: rules.rename_contigs.output
    params:
        min_identity=95,
        min_repeat=90,
        min_coverage=90,
        max_distance=20
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "fix_circular_viral_contigs_{sample}.log")
    message: "Correcting likely viral concatemers and terminal repeats in {wildcards.sample}"
    shell:
        """
        (date &&
        python ./scripts/correct_phage_contigs.py --threads {threads} --fasta {input} \
            --min_identity {params.min_identity} --min_repeat {params.min_repeat} --min_coverage {params.min_coverage} --max_distance {params.max_distance} \
            --out_fasta {output.corrected} --out_concatemer_report {output.concatemer_report} --out_dtr_report {output.dtr_report} &&
        date) &> {log}
        """

# CheckV to assess completeness
rule checkv:
    output: 
        checkv_quality = os.path.join(RESULTS_DIR, "{sample}", "checkv", "quality_summary.tsv"),
    input: 
        db = "/work/river/Databases/checkv-db-v1.5",
        assembly = rules.fix_circular_viral_contigs_per_sample.output.corrected
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_checkv.log")
    message: "Running the first CheckV per assembly"
    shell:
        """
        (
            date
            if grep -q '^>' {input.assembly:q}; then
                checkv end_to_end -t {threads} -d {input.db:q} {input.assembly:q} $(dirname {output.checkv_quality:q})
            else
                echo "No viral contigs; skipping CheckV."
                mkdir -p $(dirname {output.checkv_quality:q})
                printf 'contig_id\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\n' > {output.checkv_quality:q}
            fi
            date
        ) > {log:q} 2>&1
        """

# Produce complete per-sample and global reports, then select samples for annotation.
checkpoint viral_sample_report:
    input:
        fasta=expand(rules.fix_circular_viral_contigs_per_sample.output.corrected, sample=PHAGES_LIST),
        flye=expand(rules.assembly_reads_flye.output.info, sample=PHAGES_LIST),
        concatemer=expand(rules.fix_circular_viral_contigs_per_sample.output.concatemer_report, sample=PHAGES_LIST),
        dtr=expand(rules.fix_circular_viral_contigs_per_sample.output.dtr_report, sample=PHAGES_LIST),
        checkv=expand(rules.checkv.output.checkv_quality, sample=PHAGES_LIST),
        genomad=expand(rules.genomad.output.summary, sample=PHAGES_LIST)
    output:
        report=directory(os.path.join(RESULTS_DIR, "reports", "viral_contigs")),
        stats=expand(os.path.join(RESULTS_DIR, "{sample}", "assembly_stats.tsv"), sample=PHAGES_LIST)
    params:
        results_dir=RESULTS_DIR,
        sample_args=lambda wc, input: [value
            for entries in zip(PHAGES_LIST, input.fasta, input.flye, input.concatemer, input.dtr, input.checkv, input.genomad)
            for value in ("--sample-input", *entries)]
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "viral_sample_report.log")
    shell:
        """python ./scripts/viral_contig_report.py {params.sample_args:q} --results-dir {params.results_dir:q} --output-dir {output.report:q} > {log:q} 2>&1"""
