# running geNomad to check for viral sequences and their taxonomy
rule genomad_candidate:
    output:
        fasta=os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus.fna"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus_summary.tsv")
    input: 
        assembly = candidate_assembly,
        db = "/work/river/Databases/genomad_db/genomad_marker_metadata.tsv"
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "genomad.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{assembler}_genomad.log")
    message: "Running geNomad"
    shell:
        "(date && genomad end-to-end --threads {threads} --enable-score-calibration --composition virome --max-fdr 0.05 {input.assembly} $(dirname $(dirname {output.fasta})) $(dirname {input.db}) && date) &> {log}"

# keeping viral contigs longer than 2 kbp
rule keep_long_viral_contigs_candidate:
    output: os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "viral_above_2_kbp.fna")
    input: rules.genomad_candidate.output.fasta
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{assembler}_keep_long_viral_contigs.log")
    shell:
        """(date && seqtk seq -L 2000 {input} > {output} && date) &> {log}"""

# renaming contigs with sample name to avoid duplicates in downstream analyses
rule rename_contigs_candidate:
    output: os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "genomad", "geNomad_assembly", "assembly_summary", "viral_above_2_kbp_renamed.fna")
    input: rules.keep_long_viral_contigs_candidate.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{assembler}_rename_contigs.log")
    shell:
        """(date && awk -v sample={wildcards.sample} '/^>/ {{print ">" sample "_" substr($0, 2); next}} {{print}}' {input} > {output} && date) &> {log}"""

# autoblast to detect potential duplications of the phage (concatemers)
rule fix_circular_viral_contigs_per_sample_candidate:
    output:
        corrected=os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "circularisation", "circular_viruses.fasta"),
        concatemer_report=os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "circularisation", "breaking_concatemers_report.csv"),
        dtr_report=os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "circularisation", "breaking_dtr_report.csv")
    input: rules.rename_contigs_candidate.output
    params:
        min_identity=95,
        min_repeat=90,
        min_coverage=90,
        max_distance=20
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "fix_circular_viral_contigs_{sample}_{assembler}.log")
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
        checkv_quality = os.path.join(RESULTS_DIR, "{sample}", "{assembler}", "checkv", "quality_summary.tsv"),
    input: 
        db = "/work/river/Databases/checkv-db-v1.5",
        assembly = rules.fix_circular_viral_contigs_per_sample_candidate.output.corrected
    wildcard_constraints: assembler="flye|spades|autocycler"
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_{assembler}_checkv.log")
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

# Evaluate the primary assembly before adding any Autocycler jobs to the DAG.
checkpoint primary_viral_quality:
    input: primary_quality
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "primary_quality.tsv")
    shell: "cp {input:q} {output:q}"

# Publish one selected set of contigs for every downstream analysis.
rule select_viral_assembly:
    input: unpack(selected_candidate_inputs)
    output:
        corrected=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "circular_viruses.fasta"),
        concatemer_report=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "breaking_concatemers_report.csv"),
        dtr_report=os.path.join(RESULTS_DIR, "{sample}", "circularisation", "breaking_dtr_report.csv"),
        checkv_quality=os.path.join(RESULTS_DIR, "{sample}", "checkv", "quality_summary.tsv"),
        fasta=os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus.fna"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "genomad", "geNomad_assembly", "assembly_summary", "assembly_virus_summary.tsv"),
        assembly=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "assembly.fasta"),
        flye_info=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "flye_info.tsv"),
        assembler=os.path.join(RESULTS_DIR, "{sample}", "assembly_selection", "assembler.txt")
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
        """

# Reuse database mappings; coverage is measured on corrected viral sequences.
rule mapped_assembly_info:
    input:
        assembly=rules.select_viral_assembly.output.assembly,
        corrected=rules.select_viral_assembly.output.corrected,
        concatemer=rules.select_viral_assembly.output.concatemer_report,
        flye_info=rules.select_viral_assembly.output.flye_info,
        bam=os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam"),
        filtered_bam=os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam"),
        bai=os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam.bai"),
        filtered_bai=os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam.bai")
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_info.tsv")
    params: script=workflow.basedir + "/scripts/assembly_info.py"
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_info.log")
    shell:
        """python {params.script:q} --sample {wildcards.sample:q} --fasta {input.assembly:q} --corrected {input.corrected:q} --concatemer-report {input.concatemer:q} --bam {input.bam:q} --filtered-bam {input.filtered_bam:q} --flye-info {input.flye_info:q} --output {output:q} > {log:q} 2>&1"""

# Produce the contig table and one summary row for this sample after CheckV.
rule viral_report:
    input:
        fasta=rules.select_viral_assembly.output.corrected,
        assembly_info=rules.mapped_assembly_info.output,
        assembler=rules.select_viral_assembly.output.assembler,
        concatemer=rules.select_viral_assembly.output.concatemer_report,
        dtr=rules.select_viral_assembly.output.dtr_report,
        checkv=rules.select_viral_assembly.output.checkv_quality,
        genomad=rules.select_viral_assembly.output.summary
    output:
        contigs=os.path.join(RESULTS_DIR, "{sample}", "assembly_stats.tsv"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "viral_sample_report.tsv")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_viral_report.log")
    shell:
        """python ./scripts/viral_contig_report.py --sample {wildcards.sample:q} --assembler "$(cat {input.assembler:q})" --fasta {input.fasta:q} --assembly-info {input.assembly_info:q} --concatemer-report {input.concatemer:q} --dtr-report {input.dtr:q} --checkv {input.checkv:q} --genomad {input.genomad:q} --output {output.contigs:q} --summary-output {output.summary:q} > {log:q} 2>&1"""

# Concatenate per-sample reports before selecting samples for annotation.
checkpoint viral_report_global:
    input:
        contigs=expand(rules.viral_report.output.contigs, sample=PHAGES_LIST),
        summaries=expand(rules.viral_report.output.summary, sample=PHAGES_LIST)
    output: report=directory(os.path.join(RESULTS_DIR, "reports", "viral_contigs"))
    params:
        contig_header='Viral contig\tSample\tAssembler\tTotal bp\tTotal corrected bp\tCoverage\tCoverage without bacteria\tFlye circularity\tNumber of concatemers broken\tDTR length removed\tCheckV gene count\tCheckV viral genes\tCheckV host genes\tCheckV quality\tMIUVIG quality\tCheckV completeness\tCheckV contamination\tgeNomad provirus\tgeNomad taxonomy',
        sample_header='Sample\tViral contigs\tTotal bp\tTotal corrected bp\tStatus'
    log: os.path.join(RESULTS_DIR, "logs", "viral_report_global.log")
    shell:
        """
        (
            mkdir -p {output.report:q}
            printf '%s\n' {params.contig_header:q} > {output.report:q}/viral_contigs.tsv
            for file in {input.contigs:q}; do
                awk 'FNR > 1' "$file" >> {output.report:q}/viral_contigs.tsv
            done
            printf '%s\n' {params.sample_header:q} > {output.report:q}/samples.tsv
            for file in {input.summaries:q}; do
                awk 'FNR > 1' "$file" >> {output.report:q}/samples.tsv
            done
        ) > {log:q} 2>&1
        """
