# CheckV to assess completeness
rule checkv:
    output: 
        checkv_quality = os.path.join(RESULTS_DIR, "{sample}", "checkv", "quality_summary.tsv"),
    input: 
        db = "/work/river/Databases/checkv-db-v1.5",
        assembly = rules.fix_circular_viral_contigs_per_sample.output.corrected
    conda: os.path.join(ENV_DIR, "checkv.yaml")
    threads: 4
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_checkv.log")
    message: "Running the first CheckV per assembly"
    shell:
        """(date && checkv end_to_end -t {threads} -d {input.db:q} {input.assembly} $(dirname {output.checkv_quality}) && date) &> {log}"""

# Install into an explicit venv retained alongside the pinned source release.
rule phageterm_preparation:
    output: installation=directory(os.path.join(RESULTS_DIR, "software", "PhageTerm", "phagetermvirome-4.3"))
    log: os.path.join(RESULTS_DIR, "logs", "phageterm_preparation.log")
    conda: os.path.join(ENV_DIR, "phageterm.yaml")
    message: "Preparing PhageTermVirome"
    shell:
        """(date &&
        mkdir -p {output.installation:q} &&
        wget -O {output.installation:q}/source.tar.gz https://files.pythonhosted.org/packages/ef/89/50321c714580c79d431cd9eb12aa62dc49e6f44afbe4e3efae282c9138ff/phagetermvirome-4.3.tar.gz &&
        tar -xzf {output.installation:q}/source.tar.gz --strip-components=1 -C {output.installation:q} &&
        python -m venv {output.installation:q}/.venv &&
        export VIRTUAL_ENV="$(realpath {output.installation:q}/.venv)" &&
        export PATH="$VIRTUAL_ENV/bin:$PATH" &&
        export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True &&
        export POETRY_VIRTUALENVS_CREATE=false &&
        cd {output.installation:q} &&
        poetry install --only main --no-interaction &&
        export PYTHONPATH="$PWD/phagetermvirome" &&
        "$VIRTUAL_ENV/bin/phageterm" --help &&
        date) &> {log:q}"""


rule phageterm:
    output: os.path.join(RESULTS_DIR, "{sample}", "phageterm", "Analysis_PhageTerm_report.pdf")
    input:
        reads=rules.preprocess_reads_porechop.output,
        virus=rules.fix_circular_viral_contigs_per_sample.output.corrected,
        installation=rules.phageterm_preparation.output.installation
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_phageterm.log")
    conda: os.path.join(ENV_DIR, "phageterm.yaml")
    threads: 8
    shell:
        """(date &&
        executable="$(realpath {input.installation:q}/.venv/bin/phageterm)" &&
        export PYTHONPATH="$(realpath {input.installation:q}/phagetermvirome)" &&
        reads="$(realpath {input.reads:q})" &&
        virus="$(realpath {input.virus:q})" &&
        mkdir -p "$(dirname {output:q})" &&
        cd "$(dirname {output:q})" &&
        "$executable" -c {threads} -f "$reads" -r "$virus" -s 5 --report_title Analysis &&
        date) &> {log:q}"""
