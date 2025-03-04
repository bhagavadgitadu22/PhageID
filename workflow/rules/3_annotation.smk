rule db_pharokka:
    output: directory(os.path.join(RESULTS_DIR, "dbs", "pharokka_v1.4.0_databases"))
    log: os.path.join(RESULTS_DIR, "logs", "pharokka_db.log")
    shell:
        """(date && 
        cd $(dirname {output}) &&
        wget "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz" &&
        tar -xzf pharokka_v1.4.0_databases.tar.gz && 
        date) &> {log}"""

# Protein annotation of a phage
rule pharokka_phage:
    output: 
	    gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
	    gff = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"),
	    dnaapler = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "dnaapler", "dnaapler_reoriented.fasta"),
	    faa = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "phanotate.faa")
    input: 
        virus = rules.filtered_assembly_flye.output,
        db = rules.db_pharokka.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka.py --force -t {threads} -d {input.db} -i {input.virus} --dnaapler -o $(dirname {output.gbk}) && date) &> {log}"""

rule pharokka_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "pharokka", "plots", "{sample}_annotated_by_pharokka.png")
    input: 
        virus = rules.filtered_assembly_flye.output,
        pharokka_gbk = rules.pharokka_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && pharokka_plotter.py -i {input.virus} -n $(echo {output} | sed 's/.png//') -o $(dirname {input.pharokka_gbk}) && date) &> {log}"""

# Alternative annotation with genotate
rule genotate_calling:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate.faa")
    input: rules.pharokka_phage.output.dnaapler
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_calling.log")
    conda: os.path.join(ENV_DIR, "genotate.yaml")
    threads: 4
    shell:
        """(date && genotate.py {input} -o {output} --format faa && date) &> {log}"""

rule genotate_reformatting:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_reformatted.faa")
    input: rules.genotate_calling.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_reformatting.log")
    shell:
        """(date && sed 's/#//' {input} | sed 's/+//' > {output} && date) &> {log}"""

rule genotate_pharokka:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "pharokka_proteins_full_merged_output.tsv")
    input: 
        genes = rules.genotate_reformatting.output,
        db = rules.db_pharokka.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_annotation.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka_proteins.py --force -t {threads} -d {input.db} -i {input.genes} -o $(dirname {output}) && date) &> {log}"""

rule genotate_pharokka_gff:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "genotate_pharokka.gff")
    input: 
        metadata = rules.phage_contig_info.output,
        genes = rules.genotate_pharokka.output,
        fasta = rules.pharokka_phage.output.dnaapler
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_gff_output.log")
    shell:
        """(date && python scripts/reformat_pharokka_genotate.py {wildcards.sample} $(tail -n +2 {input.metadata} | cut -f 3) {input.genes} {input.fasta} {output} && date) &> {log}"""

rule genotate_pharokka_plot:
    output: os.path.join(RESULTS_DIR, "{sample}", "genotate", "genotate_annotation", "{sample}_annotated_by_genotate_pharokka.png")
    input:
        metadata = rules.phage_contig_info.output,
        gff_genotate = rules.genotate_pharokka_gff.output,
	    gff_pharokka = rules.pharokka_phage.output.gff
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_genotate_gff_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && python scripts/plotting_pharokka_genotate.py {wildcards.sample} $(tail -n +2 {input.metadata} | cut -f 3) {input.gff_genotate} {input.gff_pharokka} {output} && date) &> {log}"""

# DefenseFinder
rule defenseFinder_update:
    output: os.path.join(RESULTS_DIR, "logs", "defensefinder_update.txt")
    log: os.path.join(RESULTS_DIR, "logs", "defensefinder_update.log")
    conda: os.path.join(ENV_DIR, "defensefinder.yaml")
    shell:
        """(date && defense-finder update && date) > {output}"""

rule antiDefenseFinder:
    output: os.path.join(RESULTS_DIR, "{sample}", "defenseFinder", "phanotate_defense_finder_systems.tsv")
    input: 
        pharokka = rules.pharokka_phage.output.faa,
        update = rules.defenseFinder_update.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_defenseFinder.log")
    conda: os.path.join(ENV_DIR, "defensefinder.yaml")
    threads: 4    
    shell:
        """(date && 
        defense-finder run -w {threads} --preserve-raw --antidefensefinder -o $(dirname {output}) {input.pharokka} &&
        date) &> {log}"""
