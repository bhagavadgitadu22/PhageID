# preparing theBIGbam annotation file
rule reformat_empathi_results:
    output: os.path.join(RESULTS_DIR, "{sample}", "empathi", "viruses", "predictions_viruses_reformatted.csv")
    input: rules.combine_empathi_with_sublyme.output
    log: os.path.join(RESULTS_DIR, "logs", "reformat_empathi_results_{sample}.log")
    message: "Reformatting empathi results"
    shell:
        """(date &&
        python ./scripts/reformat_feature_annotations.py --input {input:q} --output {output:q} &&
        date) 2>&1 | tee {log:q}"""

rule reformat_checkamg_results:
    output: os.path.join(RESULTS_DIR, "{sample}", "checkamg", "results", "final_results_reformatted.csv")
    input:
        annotations=rules.checkamg.output,
        proteins=rules.pharokka_phage.output.faa,
        raw_proteins=rules.pharokka_phage.output.faa_raw
    log: os.path.join(RESULTS_DIR, "logs", "reformat_checkamg_results_{sample}.log")
    message: "Reformatting CheckAMG results for GFF enrichment"
    shell:
        """(date &&
        python ./scripts/reformat_feature_annotations.py --input {input.annotations:q} --output {output:q} \
            --checkamg-pharokka-fasta {input.proteins:q} --checkamg-raw-fasta {input.raw_proteins:q} &&
        date) 2>&1 | tee {log:q}"""

rule reformat_antidefence_results:
    output: os.path.join(RESULTS_DIR, "{sample}", "defenseFinder", "prodigal-gv_defense_finder_genes_reformatted.csv")
    input: rules.antiDefenseFinder.output.genes
    log: os.path.join(RESULTS_DIR, "logs", "reformat_antiDefenseFinder_{sample}.log")
    message: "Reformatting antiDefenseFinder gene results for GFF enrichment"
    shell:
        """(date &&
        python ./scripts/reformat_feature_annotations.py --input {input:q} --id-column hit_id --output {output:q} &&
        date) 2>&1 | tee {log:q}"""

# thebigbam add-contig-annotations rules
rule phold_gff:
    input: rules.phold_phage.output.gbk
    output: os.path.join(RESULTS_DIR, "{sample}", "phold", "phold.gff")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "phold_gff_{sample}.log")
    message: "Converting Phold GenBank annotations to GFF3"
    shell:
        """(date &&
        python ./scripts/genbank_to_gff.py --input {input:q} --output {output:q} &&
        date) &> {log:q}"""

rule phold_with_empathi:
    output: os.path.join(RESULTS_DIR, "{sample}", "annotations_on_viruses_enriched", "phold_with_empathi.gff")
    input:
        phold = rules.phold_gff.output,
        empathi = rules.reformat_empathi_results.output
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "phold_with_empathi_{sample}.log")
    message: "Adding Empathi predictions to Phold annotations"
    shell:
        """(date && 
        python ./scripts/enrich_annotations.py --annotations {input.phold:q} --csv {input.empathi:q} \
            --prefix empathi_ --output {output:q} &&
        date) &> {log:q}"""

rule phold_with_checkamg:
    output: os.path.join(RESULTS_DIR, "{sample}", "annotations_on_viruses_enriched", "phold_with_empathi_checkamg.gff")
    input:
        gff = rules.phold_with_empathi.output,
        checkamg = rules.reformat_checkamg_results.output
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "phold_with_checkamg_{sample}.log")
    message: "Adding CheckAMG predictions to enriched Phold annotations"
    shell:
        """(date && 
        python ./scripts/enrich_annotations.py --annotations {input.gff:q} --csv {input.checkamg:q} \
            --prefix checkamg_ --output {output:q} &&
        date) &> {log:q}"""

rule combined_viral_annotations_per_sample:
    output: os.path.join(RESULTS_DIR, "{sample}", "annotations_on_viruses_enriched", "viral_annotations_enriched.gff")
    input:
        gff=rules.phold_with_checkamg.output,
        antidefense=rules.reformat_antidefence_results.output
    conda: os.path.join(ENV_DIR, "thebigbam.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "combined_viral_annotations_{sample}.log")
    message: "Adding antiDefenseFinder predictions to enriched viral annotations"
    shell:
        """(date &&
        python ./scripts/enrich_annotations.py --annotations {input.gff:q} --csv {input.antidefense:q} \
            --prefix antidefensefinder_ --output {output:q} &&
        date) &> {log:q}"""
