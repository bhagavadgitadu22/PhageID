# protein annotation of a phage
rule pharokka_phage:
    output: 
        gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
        gff = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"),
        functions = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka_cds_functions.tsv"),
        faa = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "prodigal-gv.faa"),
        faa_raw = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "prodigal-gv_raw.faa")
    input: 
        virus = rules.select_viral_assembly.output.corrected,
        db = "/work/river/Databases/pharokka_db"
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka.py --force -t {threads} -d {input.db} -i {input.virus} -g prodigal-gv --keep_raw_prodigal -o $(dirname {output.gbk}) && date) &> {log}"""

rule pharokka_plot:
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "pharokka", "plots"))
    input: rules.pharokka_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && python ./scripts/plot_pharokka.py --genbank {input:q} --output-dir {output:q} && date) &> {log:q}"""

# Annotation with phold in addition
rule phold_phage:
    output:
        gbk = os.path.join(RESULTS_DIR, "{sample}", "phold", "phold.gbk")
    input: 
        gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
        db = "/work/river/Databases/phold_db"
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phold.log")
    conda: os.path.join(ENV_DIR, "phold.yaml")
    threads: 8
    shell:
        """(date && phold run -t {threads} --force -i {input.gbk} -d {input.db} --cpu -o $(dirname {output.gbk}) && date) &> {log}"""

# empathi on viral contigs
rule empathi_install:
    output: directory(os.path.join(RESULTS_DIR, "software", "empathi"))
    conda: os.path.join(ENV_DIR, "empathi.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "empathi_install.log")
    params:
        revision="7e9cccc00b29e11db1c0bfa066167e8026ef981c"
    message: "Installing empathi"
    shell:
        """(date && 
        mkdir -p {output:q} &&
        GIT_LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/AlexandreBoulay/empathi {output:q}/empathi &&
        git -C {output:q}/empathi lfs install --local &&
        git -C {output:q}/empathi checkout --detach {params.revision:q} &&
        git -C {output:q}/empathi lfs pull &&
        git -C {output:q}/empathi rev-parse HEAD > {output:q}/revision.txt &&
        python -m venv {output:q}/venv &&
        {output:q}/venv/bin/python -m pip install --no-cache-dir -r {output:q}/empathi/requirements.txt &&
        {output:q}/venv/bin/python -m pip install --no-cache-dir --upgrade \
            torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128 &&
        {output:q}/venv/bin/python -m pip check &&
        date) &> {log}"""

rule empathi:
    output: os.path.join(RESULTS_DIR, "{sample}", "empathi", "viruses", "predictions_viruses.csv")
    input: 
        pharokka = rules.pharokka_phage.output.faa,
        installation = rules.empathi_install.output
    log: os.path.join(RESULTS_DIR, "logs", "empathi_{sample}.log")
    conda: os.path.join(ENV_DIR, "empathi.yaml")
    params:
        hf_home=os.path.join(RESULTS_DIR, "software", "huggingface")
    threads: 4
    message: "Running empathi on all viral contigs"
    shell:
        """(date && 
        mkdir -p $(dirname $(dirname {output})) &&
        mkdir -p {params.hf_home:q}/hub &&
        export HF_HOME={params.hf_home:q} HF_HUB_CACHE={params.hf_home:q}/hub &&
        {input.installation:q}/venv/bin/python {input.installation:q}/empathi/src/empathi/empathi.py {input.pharokka:q} viruses \
            --models_folder {input.installation:q}/empathi/models --output_folder $(dirname $(dirname {output})) --threads {threads} --confidence 0.50 &&
        date) &> {log}"""

rule sublyme:
    output: os.path.join(RESULTS_DIR, "{sample}", "sublyme", "sublyme_predictions.csv")
    input: rules.pharokka_phage.output.faa
    log: os.path.join(RESULTS_DIR, "logs", "sublyme_{sample}.log")
    conda: os.path.join(ENV_DIR, "sublyme.yaml")
    threads: 4
    message: "Running sublyme on all viral contigs"
    shell:
        """(date &&
        mkdir -p $(dirname {output}) &&
        sublyme {input} --output_folder $(dirname {output}) --threads {threads} &&
        date) &> {log:q}"""

rule combine_empathi_with_sublyme:
    output: os.path.join(RESULTS_DIR, "{sample}", "empathi", "viruses", "predictions_viruses_with_sublyme.csv")
    input:
        empathi=rules.empathi.output,
        sublyme=rules.sublyme.output
    log: os.path.join(RESULTS_DIR, "logs", "combine_empathi_with_sublyme_{sample}.log")
    message: "Combining Empathi and Sublyme predictions"
    shell:
        """(date &&
        python ./scripts/combine_empathi_with_sublyme.py --empathi {input.empathi:q} --sublyme {input.sublyme:q} --output {output:q} &&
        date) > {log:q} 2>&1"""

# Report coding density and Pharokka function counts for each selected viral contig.
rule viral_annotations_report:
    input:
        fasta=rules.select_viral_assembly.output.corrected,
        gff=rules.pharokka_phage.output.gff,
        functions=rules.pharokka_phage.output.functions
    output: os.path.join(RESULTS_DIR, "{sample}", "pharokka", "viral_annotations_report.tsv")
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_viral_annotations_report.log")
    shell:
        """python ./scripts/viral_annotations_report.py --sample {wildcards.sample:q} --fasta {input.fasta:q} --gff {input.gff:q} --functions {input.functions:q} --output {output:q} > {log:q} 2>&1"""

# Concatenate reports for samples with retained viral contigs.
rule viral_annotations_report_global:
    input:
        report=lambda wc: checkpoints.viral_report_global.get().output.report,
        reports=active_sample_files(rules.viral_annotations_report.output[0])
    output: os.path.join(RESULTS_DIR, "reports", "viral_annotations", "viral_contigs.tsv")
    params:
        header="Viral_contig\tSample\tCoding_density_percent\tCDS\tconnector\tDNA, RNA and nucleotide metabolism\thead and packaging\tintegration and excision\tlysis\tmoron, auxiliary metabolic gene and host takeover\tother\ttail\ttranscription regulation\tunknown function\ttRNAs\tCRISPRs\ttmRNAs\tVFDB_Virulence_Factors\tCARD_AMR_Genes"
    log: os.path.join(RESULTS_DIR, "logs", "viral_annotations_report_global.log")
    shell:
        """(
            mkdir -p $(dirname {output:q})
            printf '%s\\n' {params.header:q} > {output:q}
            for file in {input.reports:q}; do
                awk 'FNR > 1' "$file" >> {output:q}
            done
        ) > {log:q} 2>&1"""
