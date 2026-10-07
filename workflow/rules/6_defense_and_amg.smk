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

# checkAMG on viruses for metabolic genes
rule checkamg_install:
    output:
        executable=os.path.join(RESULTS_DIR, "software", "checkamg", "venv", "bin", "checkamg")
    params:
        prefix=os.path.join(RESULTS_DIR, "software", "checkamg"),
        repository=config["checkamg"]["repository"],
        revision=config["checkamg"]["revision"]
    conda: os.path.join(ENV_DIR, "checkamg.yaml")
    log: os.path.join(RESULTS_DIR, "logs", "checkamg_install.log")
    message: "Installing pinned CheckAMG"
    shell:
        """
        (date &&
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}" &&
        python -c "import sqlite3" &&
        mkdir -p {params.prefix:q} &&
        if [ ! -e {params.prefix:q}/source ]; then
            git clone {params.repository:q} {params.prefix:q}/source
        elif [ ! -e {params.prefix:q}/source/.git ]; then
            echo "CheckAMG source directory exists but is not a Git checkout: "{params.prefix:q}/source >&2
            exit 1
        fi &&
        git -C {params.prefix:q}/source checkout --detach {params.revision:q} &&
        python -m venv {params.prefix:q}/venv &&
        {params.prefix:q}/venv/bin/python -m pip install --upgrade pip setuptools wheel packaging uv &&
        {params.prefix:q}/venv/bin/uv pip install --python {params.prefix:q}/venv/bin/python torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu &&
        {params.prefix:q}/venv/bin/uv pip install --python {params.prefix:q}/venv/bin/python torch_geometric &&
        ({params.prefix:q}/venv/bin/uv cache clean torch_scatter || true) &&
        {params.prefix:q}/venv/bin/uv pip install --python {params.prefix:q}/venv/bin/python torch_scatter --no-build-isolation --no-binary torch_scatter --refresh &&
        {params.prefix:q}/venv/bin/uv pip install --python {params.prefix:q}/venv/bin/python faiss-cpu &&
        {params.prefix:q}/venv/bin/uv pip install --python {params.prefix:q}/venv/bin/python {params.prefix:q}/source &&
        {params.prefix:q}/venv/bin/python -c "import sqlite3; import torch_scatter" &&
        {output.executable:q} --version &&
        date) &> {log:q}
        """

rule checkamg:
    input:
        pharokka = rules.pharokka.output.faa_raw,
        db = config["checkamg"]["database"],
        executable = rules.checkamg_install.output.executable
    output: os.path.join(RESULTS_DIR, "viruses", "protein_annotation", "checkamg", "results", "final_results.tsv")
    log: os.path.join(RESULTS_DIR, "logs", "checkamg_viruses.log")
    conda: os.path.join(ENV_DIR, "checkamg.yaml")
    threads: config['checkamg']['threads']    
    message: "Global checkAMG for viral contigs"
    shell:
        """(date && 
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}" &&
        export PATH="$(dirname {input.executable:q}):$PATH" &&
        {input.executable:q} annotate -t {threads} --input-type prot -p {input.pharokka:q} -d {input.db:q} -o $(dirname {output:q}) &&
        date) &> {log}"""

