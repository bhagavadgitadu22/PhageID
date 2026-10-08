# preparing theBIGbam annotation file
rule reformat_empathi_results:
    input: rules.combine_empathi_with_sublyme.output
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "empathi", "viruses", "predictions_viruses_reformatted.csv")
    log: os.path.join(RESULTS_DIR, "logs", "reformat_empathi_results.log")
    message: "Reformatting empathi results"
    params:
        converter=srcdir("../../scripts/reformat_feature_annotations.py")
    shell:
        """(date &&
        python {params.converter:q} --input {input:q} --output {output:q} &&
        date) 2>&1 | tee {log:q}"""

rule reformat_checkamg_results:
    input:
        annotations=rules.checkamg.output,
        proteins=rules.pharokka.output.faa
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "checkamg", "results", "final_results_reformatted.csv")
    log: os.path.join(RESULTS_DIR, "logs", "reformat_checkamg_results.log")
    message: "Reformatting CheckAMG results for GFF enrichment"
    params:
        converter=srcdir("../../scripts/reformat_feature_annotations.py")
    shell:
        """(date &&
        python {params.converter:q} --input {input.annotations:q} --output {output:q} \
            --checkamg-pharokka-fasta {input.proteins:q} &&
        date) 2>&1 | tee {log:q}"""

rule reformat_antidefence_results:
    # Use file-based dependency inference here: an explicit rule reference makes
    # Snakemake 7 cluster workers traverse the checkpoint merge even when its
    # outputs exist and the temporary shards have already been removed.
    input: str(rules.antiDefenseFinder_viruses.output.genes)
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "defenseFinder", "prodigal-gv_defense_finder_genes_reformatted.csv")
    log: os.path.join(RESULTS_DIR, "logs", "reformat_antidefensefinder_viruses.log")
    message: "Reformatting antiDefenseFinder gene results for GFF enrichment"
    params:
        converter=srcdir("../../scripts/reformat_feature_annotations.py")
    shell:
        """(date &&
        python {params.converter:q} --input {input:q} --id-column hit_id --output {output:q} &&
        date) 2>&1 | tee {log:q}"""

# thebigbam add-contig-annotations rules
rule pharokka_with_empathi:
    input:
        pharokka = rules.pharokka.output.gff,
        empathi = rules.reformat_empathi_results.output
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "annotations_on_viruses_enriched", "pharokka_with_empathi.gff")
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: config['thebigbam_mapping']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "pharokka_with_empathi.log")
    message: "Adding Empathi predictions to Pharokka annotations"
    shell:
        """(date && 
        thebigbam add-contig-annotations -g {input.pharokka} --csv {input.empathi} \
            --match-by feature_type,locus_tag --prefix empathi_ --keep-multiple -o {output} &&
        date) &> {log}"""

rule pharokka_with_checkamg:
    input:
        gff = rules.pharokka_with_empathi.output,
        checkamg = rules.reformat_checkamg_results.output
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "annotations_on_viruses_enriched", "pharokka_with_empathi_checkamg.gff")
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: config['thebigbam_mapping']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "pharokka_with_checkamg.log")
    message: "Adding CheckAMG predictions to enriched Pharokka annotations"
    shell:
        """(date && 
        thebigbam add-contig-annotations -g {input.gff} --csv {input.checkamg} \
            --match-by feature_type,locus_tag --prefix checkamg_ --keep-multiple -o {output} &&
        date) &> {log}"""

rule combined_viral_annotations:
    input:
        gff=rules.pharokka_with_checkamg.output,
        antidefense=rules.reformat_antidefence_results.output
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "annotations_on_viruses_enriched", "viral_annotations_enriched.gff")
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    threads: config['thebigbam_mapping']['threads']
    log: os.path.join(RESULTS_DIR, "logs", "combined_viral_annotations.log")
    message: "Adding antiDefenseFinder predictions to enriched viral annotations"
    shell:
        """(date &&
        thebigbam add-contig-annotations -g {input.gff} --csv {input.antidefense} \
            --match-by feature_type,locus_tag --prefix antidefensefinder_ --keep-multiple -o {output} &&
        date) &> {log}"""
