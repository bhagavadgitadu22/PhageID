# Termination
rule phageterm_preparation:
    output: directory(os.path.join(RESULTS_DIR, "PhageTerm", "phagetermvirome-4.3", "phagetermvirome"))
    log: os.path.join(RESULTS_DIR, "logs", "phageterm_preparation.log")
    conda: os.path.join(ENV_DIR, "phageterm.yaml")
    message: "Preparing PhageTermVirome"
    shell:
        """(date && cd $(dirname $(dirname {output})) &&
        pip install poetry &&
        export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True &&
        wget -nc https://files.pythonhosted.org/packages/ef/89/50321c714580c79d431cd9eb12aa62dc49e6f44afbe4e3efae282c9138ff/phagetermvirome-4.3.tar.gz &&
        tar -xzf phagetermvirome-4.3.tar.gz &&
        cd phagetermvirome-4.3 &&
        poetry install &&
        date) &> {log}"""

rule phageterm:
    output: os.path.join(RESULTS_DIR, "{sample}", "phageterm", "Analysis_PhageTerm_report.pdf")
    input:
#        reads = rules.concat_reads.output,
        reads = os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq"),
        virus = rules.filtered_assembly_flye.output,
        phageterm = rules.phageterm_preparation.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phageterm.log")
    conda: os.path.join(ENV_DIR, "phageterm.yaml")
    threads: 8
    shell:
        """(date && export PYTHONPATH={input.phageterm} &&
        cd $(dirname {output}) &&
        phageterm -c {threads} -f {input.reads} -r {input.virus} -s 5 --report_title Analysis && date) &> {log}"""
        
# Mapping reads on assemblies with minimap2 for COBRA
rule minimap2_mapping:
    output: 
        bam = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_mapping.bam"),
        index = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_mapping.bam.bai"),
    input:
        assembly = rules.filtered_assembly_flye.output,
#        reads = rules.concat_reads.output
        reads = os.path.join(RESULTS_DIR, "{sample}", "reads", "{sample}.fastq")
    conda: os.path.join(ENV_DIR, "phageterm_plot.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_minimap2_mapping.log")
    shell:
        """(date && \
        minimap2 -t {threads} -ax map-ont {input.assembly} {input.reads} | samtools view -@ {threads} -bS | samtools sort -@ {threads} -o {output.bam} && \
        samtools index {output.bam} {output.index} && date) &> {log}"""

rule bedgraphs_from_mapping:
    output: 
        coverage = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_coverage.bedgraph"),
        start = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_start.bedgraph"),
        end = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_end.bedgraph"),
    input: 
        bam = rules.minimap2_mapping.output.bam
    conda: os.path.join(ENV_DIR, "phageterm_plot.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_bedgraphs_from_mapping.log")
    shell:
        """(date && python scripts/plot_coverages_like_phageterm.py {input.bam} {output.coverage} {output.start} {output.end} && date) &> {log}"""

rule install_lovis4u_linux:
    output: os.path.join(RESULTS_DIR, "logs", "lovis4u_linux.touch")
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    shell:
        """lovis4u --linux && lovis4u --get-hmms && touch {output}"""

rule phage_term_plot:
    output:
        os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "lovis4u.pdf"),
    input:
        linux = rules.install_lovis4u_linux.output,
        gbk = rules.pharokka_phage.output.gbk,
        coverage = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_coverage.bedgraph"),
        start = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_start.bedgraph"),
        end = os.path.join(RESULTS_DIR, "{sample}", "phageterm_plot", "{sample}_end.bedgraph"),
    conda: os.path.join(ENV_DIR, "lovis4u.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phageterm_plot.log")
    shell:
        """(date && lovis4u -gb {input.gbk} --set-category-colour -o $(dirname {output}) \
        -scc -gc -gc_skew -bg {input.coverage} {input.start} {input.end} -bgl "coverage" "start_reads" "end_reads" && date) &> {log}"""
