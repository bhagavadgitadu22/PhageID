# protein annotation of a phage
rule pharokka_phage:
    output: 
        gbk = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gbk"),
        gff = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "pharokka.gff"),
        faa = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "phanotate.faa"),
        faa_raw = os.path.join(RESULTS_DIR, "{sample}", "pharokka", "phanotate_raw.faa")
    input: 
        virus = rules.fix_circular_viral_contigs_per_sample.output.corrected,
        db = "/work/river/Databases/pharokka_db"
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    threads: 8
    shell:
        """(date && pharokka.py --force -t {threads} -d {input.db} -i {input.virus} -o $(dirname {output.gbk}) && date) &> {log}"""

rule pharokka_plot:
    output: directory(os.path.join(RESULTS_DIR, "{sample}", "pharokka", "plots"))
    input: rules.pharokka_phage.output.gbk
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_pharokka_plot.log")
    conda: os.path.join(ENV_DIR, "pharokka.yaml")
    shell:
        """(date && pharokka multiplot -g {input} -o {output} && date) &> {log}"""

# Annotation with phold in addition
rule phold_phage:
    output: 
        gff = os.path.join(RESULTS_DIR, "{sample}", "phold", "phold.gff"),
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
    input: 
        pharokka = rules.pharokka_phage.output.faa,
        installation = rules.empathi_install.output
    output: os.path.join(RESULTS_DIR, "{sample}", "empathi", "viruses", "predictions_viruses.csv")
    log: os.path.join(RESULTS_DIR, "logs", "empathi_{sample}.log")
    conda: os.path.join(ENV_DIR, "empathi.yaml")
    threads: 4
    message: "Running empathi on all viral contigs"
    shell:
        """(date && 
        mkdir -p $(dirname $(dirname {output})) &&
        {input.installation:q}/venv/bin/python {input.installation:q}/empathi/src/empathi/empathi.py {input.pharokka:q} viruses \
            --models_folder {input.installation:q}/empathi/models --output_folder $(dirname $(dirname {output})) --threads {threads} --confidence 0.50 &&
        date) &> {log}"""

rule sublyme:
    input: rules.pharokka_phage.output.faa
    output: os.path.join(RESULTS_DIR, "{sample}", "sublyme", "sublyme_predictions.csv")
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
    input:
        empathi=rules.empathi.output,
        sublyme=rules.sublyme.output
    output: os.path.join(RESULTS_DIR, "{sample}", "empathi", "viruses", "predictions_viruses_with_sublyme.csv")
    log: os.path.join(RESULTS_DIR, "logs", "combine_empathi_with_sublyme_{sample}.log")
    message: "Combining Empathi and Sublyme predictions"
    shell:
        """(date &&
        python ./scripts/combine_empathi_with_sublyme.py --empathi {input.empathi:q} --sublyme {input.sublyme:q} --output {output:q} &&
        date) > {log:q} 2>&1"""

