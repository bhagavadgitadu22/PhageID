# Removing remaining adapters and barcodes
rule preprocess_reads_porechop:
    output: os.path.join(RESULTS_DIR, "{sample}", "reads", "porechop.{sample}.fastq")
    input: os.path.join(READS_DIR, "{sample}")
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_preprocess_reads_porechop.log")
    threads: 10
    shell:
        """(date && porechop -t {threads} -i {input} -o {output} --discard_middle --check_reads 10000 --end_size 200 && date) &> {log}"""

# Assembling long-reads with flye
# meta recommended for extrachromosomal elements like phages or plasmids
rule assembly_reads_flye:
    output:
        assembly = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly.fasta"),
        graph = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_graph.gfa"),
        info = os.path.join(RESULTS_DIR, "{sample}", "flye", "assembly_info.txt")
    input: rules.preprocess_reads_porechop.output
    conda: os.path.join(ENV_DIR, "assembly_flye.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_reads_flye.log")
    threads: 20
    shell:
        """(date && flye -t {threads} --meta --nano-raw {input} -o $(dirname {output.assembly}) && date) &> {log}"""

# These 2 rules try to keep only complete phage contig from the assembly
# Phage contigs where picked manually when this failed (looking for long circular high coverage contigs)
rule depths_per_contig_flye:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "contigs_to_keep.txt")
    input: rules.assembly_reads_flye.output.info
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_depths_per_contig_flye.log")
    shell:
        """(date &&
        # Count the number of contigs
        contig_count=$(awk 'NR>1 {{print $1}}' {input} | wc -l)

        # Extract the sequence names of contigs that circularize
        circular_contigs=$(awk 'NR>1 && $4 == "Y" {{print $1}}' {input})
        circular_count=$(echo "$circular_contigs" | wc -l)

        # Check if there is only one contig or only one circular contig
        if [ "$contig_count" -eq 1 ]; then
            # Extract the sequence name of the single contig
            awk 'NR>1 {{print $1}}' {input} > {output}
        elif [ "$circular_count" -eq 1 ]; then
            # Write the number of the single circular contig
            echo "$circular_contigs" > {output} 
        fi &&
        date) &> {log}"""

# For a few recalcitrant phages, TruSeq DNA was produced in addition to the Nanopore reads
# Hybrid assembly was performed with unicycler and the resulting assembly replaced the output of the next rule
rule filtered_assembly_flye:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "filtered_assembly.fasta")
    input:
        assembly = rules.assembly_reads_flye.output.assembly,
        list_contigs = rules.depths_per_contig_flye.output
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_filtered_assembly_flye.log")
    shell:
        """(date && seqtk subseq {input.assembly} {input.list_contigs} > {output} &&
        sed -i "s/>/>{wildcards.sample}_/" {output} && date) &> {log}"""

# Autoblast to detect potential duplications of the phage (concatemers)
rule correct_phage_with_autoblast:
    output: 
        corrected = os.path.join(RESULTS_DIR, "{sample}", "flye", "autoblast_corrected.fasta"),
        report = os.path.join(RESULTS_DIR, "{sample}", "flye", "autoblast_report.txt")
    input: rules.filtered_assembly_flye.output
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_autoblast.log")
    shell:
        r"""
        (date && 
        # Run BLAST against itself
        blastn -query {input} -subject {input} -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen" -evalue 1e-10 > {input}.blast.tsv 2>> {log}
        # Analyze results for duplications
        python scripts/analyse_duplications.py --blast {input}.blast.tsv --fasta {input} --out_report {output.report} --out_fasta {output.corrected} &&
        date) &> {log}
        """

# Assessing final contigs
rule quality_assembly:
    output: os.path.join(RESULTS_DIR, "{sample}", "flye", "quast", "report.html")
    input: rules.correct_phage_with_autoblast.output.corrected
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_quality_assembly.log")
    shell:
        """(date && quast.py {input} -o $(dirname {output}) && date) &> {log}"""

rule phage_contig_info:
    output: os.path.join(RESULTS_DIR, "{sample}", "assembly_stats.tsv")
    input: rules.correct_phage_with_autoblast.output.corrected
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_assembly_stats.log")
    shell:
        """(date && ./scripts/contig_info.sh -m 1000 -t {input} > {output} && date) &> {log}"""
