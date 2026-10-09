rule genemap_for_vcontact3:
    output:
        proteins=os.path.join(RESULTS_DIR, "vcontact3_inputs", "proteins.faa"),
        gene_map=os.path.join(RESULTS_DIR, "vcontact3_inputs", "gene2genome.tsv"),
        lengths=os.path.join(RESULTS_DIR, "vcontact3_inputs", "genome_lengths.tsv")
    input:
        representatives=rules.viruses_dereplicated.output.fna_viruses_derep,
        genbanks=active_sample_files(os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"))
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "genemap_for_vcontact3.log")
    shell:
        """python ./scripts/prepare_vcontact3.py --representatives {input.representatives:q} --genbanks {input.genbanks:q} \
            --proteins {output.proteins:q} --gene-map {output.gene_map:q} --genome-lengths {output.lengths:q} > {log:q} 2>&1"""

rule vcontact3:
    output: directory(os.path.join(RESULTS_DIR, "vcontact3_results"))
    input:
        faa=rules.genemap_for_vcontact3.output.proteins,
        gene_map=rules.genemap_for_vcontact3.output.gene_map,
        lengths=rules.genemap_for_vcontact3.output.lengths,
        db="/work/river/Databases/vContact3_db"
    conda: os.path.join(ENV_DIR, "vcontact3.yaml")
    threads: 8
    log: os.path.join(RESULTS_DIR, "logs", "vcontact3.log")
    message: "Running vContact3 on dereplicated viruses"
    shell:
        """(date &&
        vcontact3 run --threads {threads} --db-path {input.db:q} --db-domain prokaryotes \
            --proteins {input.faa:q} --gene2genome {input.gene_map:q} --len-nucleotide {input.lengths:q} \
            --output {output:q} --exports cytoscape profiles completeness &&
        date) > {log:q} 2>&1"""
