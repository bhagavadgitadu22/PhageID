# running geNomad to check for viral sequences and their taxonomy
# readset can be filtered|all (filtered bacterial reads or all reads)
# assembler can be flye|spades|autocycler
rule genomad_candidate:
    output:
        fasta=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus.fna"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus_summary.tsv")
    input: 
        assembly = candidate_assembly,
        db = "/work/river/Databases/genomad_db/genomad_marker_metadata.tsv"
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_{assembler}_genomad.log")
    message: "Running geNomad"
    shell:
        "(date && genomad end-to-end --threads {threads} --enable-score-calibration --composition virome --max-fdr 0.05 {input.assembly} $(dirname $(dirname {output.fasta})) $(dirname {input.db}) && date) &> {log}"

# keeping viral contigs longer than 2 kbp
rule keep_long_viral_contigs_candidate:
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "viral_above_2_kbp.fna")
    input: rules.genomad_candidate.output.fasta
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_{assembler}_keep_long_viral_contigs.log")
    shell:
        """(date && seqtk seq -L 2000 {input} > {output} && date) &> {log}"""

# renaming contigs with sample name to avoid duplicates in downstream analyses
rule rename_contigs_candidate:
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "viral_above_2_kbp_renamed.fna")
    input: rules.keep_long_viral_contigs_candidate.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_{assembler}_rename_contigs.log")
    shell:
        """(date && awk -v sample={wildcards.sample} '/^>/ {{print ">" sample "_" substr($0, 2); next}} {{print}}' {input} > {output} && date) &> {log}"""

# autoblast to detect potential duplications of the phage (concatemers)
rule fix_circular_viral_contigs_per_sample_candidate:
    output:
        corrected=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "circularisation", "circular_viruses.fasta"),
        concatemer_report=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "circularisation", "breaking_concatemers_report.csv"),
        dtr_report=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "circularisation", "breaking_dtr_report.csv")
    input: rules.rename_contigs_candidate.output
    params:
        min_identity=95,
        min_repeat=90,
        min_coverage=90,
        max_distance=20
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "fix_circular_viral_contigs_{sample}_{readset}_{assembler}.log")
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
rule checkv_candidate:
    output: 
        checkv_quality = os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "checkv", "quality_summary.tsv"),
    input: 
        db = "/work/river/Databases/checkv-db-v1.5",
        assembly = rules.fix_circular_viral_contigs_per_sample_candidate.output.corrected
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{readset}_{assembler}_checkv.log")
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
                printf 'contig_id\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcontamination\n' > {output.checkv_quality:q}
            fi
            date
        ) > {log:q} 2>&1
        """

# Each quality checkpoint lets selection inspect a completed candidate.
checkpoint candidate_viral_quality:
    input: rules.checkv_candidate.output.checkv_quality
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "validated_quality.tsv")
    wildcard_constraints: assembler="flye|spades|autocycler"
    shell: "cp {input:q} {output:q}"

# Publish one selected set of contigs for every downstream analysis.
rule select_viral_assembly:
    input:
        corrected=lambda wc: selected_candidate_inputs(wc)["corrected"],
        concatemer_report=lambda wc: selected_candidate_inputs(wc)["concatemer_report"],
        dtr_report=lambda wc: selected_candidate_inputs(wc)["dtr_report"],
        checkv_quality=lambda wc: selected_candidate_inputs(wc)["checkv_quality"],
        fasta=lambda wc: selected_candidate_inputs(wc)["fasta"],
        summary=lambda wc: selected_candidate_inputs(wc)["summary"],
        assembly=lambda wc: selected_candidate_inputs(wc)["assembly"],
        used_reads=lambda wc: selected_candidate_inputs(wc)["used_reads"],
        filter_report=lambda wc: selected_candidate_inputs(wc)["filter_report"],
        flye_info=lambda wc: [selected_candidate_inputs(wc)["flye_info"]] if "flye_info" in selected_candidate_inputs(wc) else []
    output:
        corrected=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"),
        concatemer_report=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "breaking_concatemers_report.csv"),
        dtr_report=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "breaking_dtr_report.csv"),
        checkv_quality=os.path.join(RESULTS_DIR, "{sample}", "checkv", "quality_summary.tsv"),
        fasta=os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus.fna"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus_summary.tsv"),
        assembly=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "assembly.fasta"),
        flye_info=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "flye_info.tsv"),
        assembler=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "assembler.txt"),
        read_stats=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "read_stats.tsv")
    shell:
        """
        cp {input.corrected:q} {output.corrected:q}
        cp {input.concatemer_report:q} {output.concatemer_report:q}
        cp {input.dtr_report:q} {output.dtr_report:q}
        cp {input.checkv_quality:q} {output.checkv_quality:q}
        cp {input.fasta:q} {output.fasta:q}
        cp {input.summary:q} {output.summary:q}
        cp {input.assembly:q} {output.assembly:q}
        if [ "$(basename "$(dirname "$(dirname {input.corrected:q})")")" = flye ]; then
            cp "$(dirname {input.assembly:q})/assembly_info.txt" {output.flye_info:q}
        else
            printf '#seq_name\tlength\tcov.\tcirc.\n' > {output.flye_info:q}
        fi
        basename "$(dirname "$(dirname {input.corrected:q})")" > {output.assembler:q}
        readset=$(basename "$(dirname "$(dirname {input.assembly:q})")")
        python ./scripts/assembly_read_stats.py selected --filter-report {input.filter_report:q} --used-reads {input.used_reads:q} --readset "$readset" --output {output.read_stats:q}
        """
