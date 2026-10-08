# protein annotation of a phage
rule pharokka_phage:
    output: 
	    gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
	    gff = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"),
	    dnaapler = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "dnaapler", "dnaapler_reoriented.fasta"),
	    faa = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "phanotate.faa")
    input: 
        virus = rules.filtered_assembly_flye.output,
        db = config["pharokka"]["database"]
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

# Annotation with phold in addition
rule phold_phage:
    output: 
	    gbk = os.path.join(RESULTS_DIR, "{sample}", "phold", "phold.gbk")
    input: 
        gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
        db = config["phold"]["database"]
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phold.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && phold run -t {threads} --force -i {input.gbk} -d {input.db} --cpu -o $(dirname {output.gbk}) && date) &> {log}"""

rule phold_plot:
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "phold", "plots"))
    input: 
        phold_gbk = rules.phold_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phold_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && phold plot --force -i {input.phold_gbk} -o {output} && date) &> {log}"""

# empathi on viral contigs
rule empathi_install:
    output: directory(os.path.join(RESULTS_DIR, "software", "empathi"))
    conda: os.path.join(ENV_DIR, "empathi.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "empathi_install.log")
    params:
        revision=config["empathi"]["revision"]
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

# GPU
rule empathi:
    input: 
        pharokka = rules.pharokka.output.faa,
        installation = rules.empathi_install.output
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "empathi", "viruses", "predictions_viruses.csv")
    params:
        outdir = os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "empathi")
    log: os.path.join(RESULTS_DIR, "logs", "empathi_viruses.log")
    conda: os.path.join(ENV_DIR, "empathi.yaml")
    threads: config['empathi']['threads']
    message: "Running empathi on all viral contigs"
    shell:
        """(date && 
        mkdir -p {params.outdir:q} &&
        {input.installation:q}/venv/bin/python {input.installation:q}/empathi/src/empathi/empathi.py {input.pharokka:q} viruses \
            --models_folder {input.installation:q}/empathi/models --output_folder {params.outdir:q} --threads {threads} --confidence 0.50 &&
        date) &> {log}"""

rule sublyme:
    input:
        pharokka = rules.pharokka.output.faa
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "sublyme", "sublyme_predictions.csv")
    params:
        outdir = os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "sublyme")
    log: os.path.join(RESULTS_DIR, "logs", "sublyme_viruses.log")
    conda: os.path.join(ENV_DIR, "sublyme.yaml")
    threads: config['empathi']['threads']
    message: "Running sublyme on all viral contigs"
    shell:
        """(date &&
        mkdir -p {params.outdir:q} &&
        sublyme {input.pharokka:q} --output_folder {params.outdir:q} --threads {threads} &&
        date) &> {log:q}"""

rule combine_empathi_with_sublyme:
    input:
        empathi=rules.empathi.output,
        sublyme=rules.sublyme.output
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "empathi", "viruses", "predictions_viruses_with_sublyme.csv")
    log: os.path.join(RESULTS_DIR, "logs", "combine_empathi_with_sublyme.log")
    params:
        converter=srcdir("../../scripts/combine_empathi_with_sublyme.py")
    message: "Combining Empathi and Sublyme predictions"
    shell:
        """(date &&
        python {params.converter:q} --empathi {input.empathi:q} --sublyme {input.sublyme:q} --output {output:q} &&
        date) > {log:q} 2>&1"""
