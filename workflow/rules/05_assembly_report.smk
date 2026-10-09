# theBIGbam for all phages individually
rule thebigbam_mapping:
    output: 
        bam = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam"),
        bai = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}.bam.bai")
    input:
        assembly = rules.select_viral_assembly.output.corrected,
        read1 = sample_reads
    params: mapper=sample_mapper
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_mapping.log")
    message: "Running minimap2 to calculate coverage"
    shell:
        """(date &&
        if grep -q '^>' {input.assembly:q}; then
            thebigbam mapping-per-sample -t {threads} -r1 {input.read1:q} -a {input.assembly:q} --mapper {params.mapper:q} --circular -o {output.bam:q}
        else
            printf '@HD\tVN:1.6\tSO:coordinate\n@CO\ttheBIGbam:circular=true\n' | samtools view -b -o {output.bam:q} -
            samtools index {output.bam:q}
        fi &&
        date) &> {log}"""

rule thebigbam_mapping_without_bacteria:
    output: 
        bam = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam"),
        bai = os.path.join(RESULTS_DIR, "minimap2", "thebigbam", "{sample}_without_bacteria.bam.bai")
    input:
        assembly = rules.select_viral_assembly.output.corrected,
        read1 = rules.remove_bacterial_contamination.output
    params: mapper=sample_mapper
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_mapping_without_bacteria.log")
    message: "Running minimap2 to calculate coverage"
    shell:
        """(date &&
        if grep -q '^>' {input.assembly:q}; then
            thebigbam mapping-per-sample -t {threads} -r1 {input.read1:q} -a {input.assembly:q} --mapper {params.mapper:q} --circular -o {output.bam:q}
        else
            printf '@HD\tVN:1.6\tSO:coordinate\n@CO\ttheBIGbam:circular=true\n' | samtools view -b -o {output.bam:q} -
            samtools index {output.bam:q}
        fi &&
        date) &> {log}"""

# Reuse database mappings; coverage is measured on corrected viral sequences.
rule mapped_assembly_info:
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_info.tsv")
    input:
        assembly=rules.select_viral_assembly.output.assembly,
        corrected=rules.select_viral_assembly.output.corrected,
        concatemer=rules.select_viral_assembly.output.concatemer_report,
        flye_info=rules.select_viral_assembly.output.flye_info,
        bam=rule.thebigbam_mapping.output.bam,
        filtered_bam=rule.thebigbam_mapping_without_bacteria.output.bam,
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_info.log")
    shell:
        """python ./scripts/assembly_info.py --sample {wildcards.sample:q} --fasta {input.assembly:q} --corrected {input.corrected:q} --concatemer-report {input.concatemer:q} --bam {input.bam:q} --filtered-bam {input.filtered_bam:q} --flye-info {input.flye_info:q} --output {output:q} > {log:q} 2>&1"""

# Produce the contig table and one summary row for this sample after CheckV.
rule viral_report:
    input:
        fasta=rules.select_viral_assembly.output.corrected,
        assembly_info=rules.mapped_assembly_info.output,
        assembler=rules.select_viral_assembly.output.assembler,
        read_stats=rules.select_viral_assembly.output.read_stats,
        concatemer=rules.select_viral_assembly.output.concatemer_report,
        dtr=rules.select_viral_assembly.output.dtr_report,
        checkv=rules.select_viral_assembly.output.checkv_quality,
        genomad=rules.select_viral_assembly.output.summary
    output:
        contigs=os.path.join(RESULTS_DIR, "{sample}", "assembly_stats.tsv"),
        summary=os.path.join(RESULTS_DIR, "{sample}", "viral_sample_report.tsv")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_viral_report.log")
    shell:
        """python ./scripts/viral_contig_report.py --sample {wildcards.sample:q} --assembler "$(cat {input.assembler:q})" --read-stats {input.read_stats:q} --fasta {input.fasta:q} --assembly-info {input.assembly_info:q} --concatemer-report {input.concatemer:q} --dtr-report {input.dtr:q} --checkv {input.checkv:q} --genomad {input.genomad:q} --output {output.contigs:q} --summary-output {output.summary:q} > {log:q} 2>&1"""

# Concatenate per-sample reports before selecting samples for annotation.
checkpoint viral_report_global:
    input:
        contigs=expand(rules.viral_report.output.contigs, sample=PHAGES_LIST),
        summaries=expand(rules.viral_report.output.summary, sample=PHAGES_LIST)
    output: report=directory(os.path.join(RESULTS_DIR, "reports", "viral_contigs"))
    params:
        contig_header='Viral contig\tSample\tAssembler\tAssembly reads\tPercentage bacterial reads\tAll reads number\tReads used for assembly number\tTotal bp\tTotal corrected bp\tCoverage\tCoverage without bacteria\tFlye circularity\tNumber of concatemers broken\tDTR length removed\tCheckV gene count\tCheckV viral genes\tCheckV host genes\tCheckV quality\tMIUVIG quality\tCheckV completeness\tCheckV contamination\tgeNomad provirus\tgeNomad taxonomy',
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
