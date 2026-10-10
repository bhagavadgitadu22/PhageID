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
        """
        (
            if grep -q '^>' {input.assembly:q}; then
                genomad end-to-end --threads {threads} --enable-score-calibration --composition virome --max-fdr 0.05 {input.assembly:q} $(dirname $(dirname {output.fasta:q})) $(dirname {input.db:q})
            else
                echo "Empty or failed assembly; skipping geNomad."
                mkdir -p $(dirname {output.fasta:q})
                : > {output.fasta:q}
                printf 'seq_name\ttopology\ttaxonomy\n' > {output.summary:q}
            fi
        ) > {log:q} 2>&1
        """

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
            # CheckV otherwise reuses cached genes from a previous input FASTA.
            checkv_dir=$(dirname {output.checkv_quality:q})
            rm -rf -- "$checkv_dir"
            mkdir -p "$checkv_dir"
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
    input:
        quality=rules.checkv_candidate.output.checkv_quality,
        status=lambda wc: os.path.join(os.path.dirname(candidate_assembly(wc)), "assembly_status.txt")
    output:
        quality=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "validated_quality.tsv"),
        status=os.path.join(RESULTS_DIR, "{sample}", "assembly_attempts", "{readset}", "{assembler}", "validated_assembly_status.txt")
    wildcard_constraints: assembler="flye|spades|autocycler"
    shell: "cp {input.quality:q} {output.quality:q}; cp {input.status:q} {output.status:q}"

# Resolve the choice once and persist validated source paths before publication.
rule assembly_selection_manifest:
    input: sources=selection_manifest_inputs
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "selected_sources.json")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_selection_manifest.log")
    shell:
        """python ./scripts/assembly_selection_manifest.py create --sample {wildcards.sample:q} --inputs {input.sources:q} --output {output:q} > {log:q} 2>&1"""

# Static dependency: forcing publication cannot substitute checkpoint placeholders.
rule select_viral_assembly:
    input: manifest=rules.assembly_selection_manifest.output
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
        status=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "assembly_status.txt"),
        read_stats=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "read_stats.tsv")
    params:
        destinations=lambda wc, output: [value for name, path in output.items()
                                         for value in ("--destination", name, str(path))]
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_select_viral_assembly.log")
    shell:
        """python ./scripts/assembly_selection_manifest.py publish --sample {wildcards.sample:q} --manifest {input.manifest:q} {params.destinations:q} > {log:q} 2>&1"""
