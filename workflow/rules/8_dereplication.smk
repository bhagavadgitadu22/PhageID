# Copy all fasta and gff from phages of a host in a common directory
rule copy_files_genus:
    output: 
        fasta = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "flye_assemblies")),
        gff_pharokka = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "pharokka_annotations")),
        gff_genotate = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "genotate_annotations")),
        pharokka_plot = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "pharokka_plots")),
        genotate_plot = directory(os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "genotate_plots")),
    input: 
        # Match all directories starting with genus
        lambda wildcards: sorted(
            [os.path.join(RESULTS_DIR, d, "flye", "filtered_assembly.fasta") 
            for d in os.listdir(RESULTS_DIR) 
            if d.startswith(wildcards.genus + "_") 
            and os.path.exists(os.path.join(RESULTS_DIR, d, "flye", "filtered_assembly.fasta"))]
        )
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_copy_files.log")
    shell:
        """(date && mkdir -p $(dirname {output.fasta} | xargs dirname) && mkdir -p $(dirname {output.fasta}) &&
            mkdir -p {output.fasta} && mkdir -p {output.gff_pharokka} && mkdir -p {output.gff_genotate} && 
            mkdir -p {output.pharokka_plot} && mkdir -p {output.genotate_plot}
            for f in {input}; do 
                name=$(dirname $f | xargs dirname | xargs basename);
                ln -s $f {output.fasta}/${{name}}.fasta;
                
                gff_pharokka=$(echo $f | sed 's/flye\/filtered_assembly.fasta/pharokka\/pharokka.gff/');
                gff_genotate=$(echo $f | sed 's/flye\/filtered_assembly.fasta/genotate\/genotate_annotation\/genotate_pharokka.gff/');
                plot_pharokka=$(echo $f | sed "s/flye\/filtered_assembly.fasta/pharokka\/plots\/${{name}}_annotated_by_pharokka.png/");
                plot_genotate=$(echo $f | sed "s/flye\/filtered_assembly.fasta/genotate\/genotate_annotation\/${{name}}_annotated_by_genotate_pharokka.png/");
                echo $plot_pharokka
                
                ln -s $gff_pharokka {output.gff_pharokka}/pharokka_${{name}}.gff;
                ln -s $gff_genotate {output.gff_genotate}/genotate_${{name}}.gff;
                cp $plot_pharokka {output.pharokka_plot}/${{name}}_annotated_by_pharokka.png;
                cp $plot_genotate {output.genotate_plot}/${{name}}_annotated_by_genotate.png;
            done && date) &> {log}"""

# Comparison of phages' sequences relative to this host
rule fast_ani:
    output: 
        list = os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "list_sequences.txt"),
        fastani = os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "fastani.out")
    input: rules.copy_files_genus.output.fasta
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_fastani.log")
    shell:
        """(date && find {input} -name "*.fasta" > {output.list} &&
        fastANI --rl {output.list} --ql {output.list} --matrix -o {output.fastani} && date) &> {log}"""

# Comparison of phages' annotations relative to this host
rule install_lovis4u_linux:
    output: os.path.join(RESULTS_DIR, "logs", "lovis4u_linux.touch")
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    shell:
        """lovis4u --linux && lovis4u --get-hmms && touch {output}"""

rule lovis4u_pharokka:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "lovis4u_pharokka", "lovis4u.pdf")
    input: 
        linux = rules.install_lovis4u_linux.output,
        gff = rules.copy_files_genus.output.gff_pharokka
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_lovis4u_pharokka.log")
    shell:
        """(date && lovis4u -gff {input.gff} --reorient_loci --use-filename-as-id --homology-links --set-category-colour --run-hmmscan -o $(dirname {output}) && date) &> {log}"""

rule lovis4u_genotate:
    output: os.path.join(RESULTS_DIR, "analysis_per_host", "{genus}", "lovis4u_genotate", "lovis4u.pdf")
    input:
        linux = rules.install_lovis4u_linux.output,
        gff = rules.copy_files_genus.output.gff_genotate
    conda: os.path.join(ENV_DIR, "dereplication.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "{genus}_lovis4u_genotate.log")
    shell:
        """(date && lovis4u -gff {input.gff} --reorient_loci --use-filename-as-id --homology-links --set-category-colour --run-hmmscan -o $(dirname {output}) && date) &> {log}"""

# theBIGbam
rule thebigbam_mapping:
    input:
        assembly = choose_input_file,
        read1 = os.path.join(RESULTS_DIR, "cutadapt", "{assembly}_R1_cutadapt.fastq.gz"),
        read2 = os.path.join(RESULTS_DIR, "cutadapt", "{assembly}_R2_cutadapt.fastq.gz")
    output: 
        bam = os.path.join(RESULTS_DIR, "{kingdom}", "minimap2", "thebigbam", "{assembly}.bam"),
        bai = os.path.join(RESULTS_DIR, "{kingdom}", "minimap2", "thebigbam", "{assembly}.bam.bai")
    params:
        reference_bases_dealt_at_once = "16G",
        min_read_identity = config['coverm']['min_read_identity'],
        min_read_coverage = config['coverm']['min_read_coverage'],
        circular = lambda wildcards: "--circular" if "viruses" in wildcards.kingdom else ""
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: config['minimap2']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "{assembly}_mapping_for_coverage_{kingdom}.log")
    message: "Running minimap2 to calculate coverage"
    shell:
        """(date &&
        thebigbam mapping-per-sample -t {threads} \
            -r1 {input.read1:q} -r2 {input.read2:q} -a {input.assembly:q} \
            --min-read-percent-identity {params.min_read_identity} \
            --min-read-aligned-percent {params.min_read_coverage} \
            {params.circular} -o {output.bam:q} &&
        date) &> {log}"""

rule thebigbam_calculate:
    input:
        annotation = lambda wildcards: rules.combined_viral_annotations.output if "viruses" in wildcards.kingdom else rules.combine_mag_annotations.output,
        bam = expand(os.path.join(RESULTS_DIR, "{{kingdom}}", "minimap2", "thebigbam", "{assembly}.bam"), assembly=SAMPLES),
        bai = expand(os.path.join(RESULTS_DIR, "{{kingdom}}", "minimap2", "thebigbam", "{assembly}.bam.bai"), assembly=SAMPLES)
    output:
        db = os.path.join(RESULTS_DIR, "{kingdom}", "thebigbam", "thebigbam_{kingdom}.db")
    params:
        bam_dir = os.path.join(RESULTS_DIR, "{kingdom}", "minimap2", "thebigbam"),
        view = lambda wildcards: "contig" if "viruses" in wildcards.kingdom else "mag",
        min_aligned_fraction = lambda wildcards: "" if "viruses" in wildcards.kingdom else "--min_aligned_fraction 10"
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: config['thebigbam_calculate']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "thebigbam_calculate_for_{kingdom}.log")
    shell:
        """(date &&
        thebigbam calculate -t {threads} -b {params.bam_dir:q} -g {input.annotation:q} -o {output.db:q} --time --blast {params.min_aligned_fraction} --view {params.view} &&
        date) &> {log}"""
