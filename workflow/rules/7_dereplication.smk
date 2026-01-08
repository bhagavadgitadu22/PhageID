# Dereplication of the viruses
rule combine_all_viruses:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "all_viruses.fna")
    input: expand(os.path.join(RESULTS_DIR, "{sample}", "flye", "autoblast_corrected.fasta"), sample = PHAGES_LIST)
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "combine_all_viruses.log")
    message: "Combining all viruses from all samples"
    shell:
        """
        (date && mkdir -p $(dirname {output}) &&
        cat {input} > {output} && date) &> {log}
        """

rule blast_before_dereplication:
    output: os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "blastn_all_viruses.tsv")
    input: rules.combine_all_viruses.output
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "blastn_all_viruses.log")
    message: "BLAST all-vs-all of the viruses"
    shell:
        """
        (date && blastn -num_threads {threads} -query {input} -subject {input} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out {output} &&
        date) &> {log}
        """

rule ani_for_dereplication:
    output: 
        ani_results = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "ani_all_viruses.tsv"),
        clustering_results = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "clusters_all_viruses.tsv"),
    input:
        blast_results = rules.blast_before_dereplication.output, 
        fna_viruses = rules.combine_all_viruses.output
    conda: os.path.join(ENV_DIR, "viral_detection.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "ani_all_viruses.log")
    message: "ANI for dereplication of the viral dataset"
    shell:
        """
        (date && python scripts/ani_calc.py -i {input.blast_results} -o {output.ani_results} &&
        python scripts/ani_clust.py --fna {input.fna_viruses} --ani {output.ani_results} --out {output.clustering_results} --min_ani 95 --min_tcov 85 --min_qcov 0 && date) &> {log}
        """

rule viruses_dereplicated:
    output:
        list_viruses_derep = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "all_viruses_dereplicated.txt"),
        fna_viruses_derep = os.path.join(RESULTS_DIR, "combined_viruses", "dereplication", "all_viruses_dereplicated.fna"),
    input:
        fna_viruses = rules.combine_all_viruses.output,
        clustering_results = rules.ani_for_dereplication.output.clustering_results,
    conda: os.path.join(ENV_DIR, "preprocessing.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "dereplication_all_viruses.log")
    message: "Dereplication of the viral dataset"
    shell:
        """
        (date && cut -f 1 {input.clustering_results} > {output.list_viruses_derep} &&
        seqtk subseq {input.fna_viruses} {output.list_viruses_derep} > {output.fna_viruses_derep} && date) &> {log}
        """

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
        gff = expand(os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"), sample = PHAGES_LIST)
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
        gff = expand(os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "genotate_pharokka.gff"), sample = PHAGES_LIST)
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
