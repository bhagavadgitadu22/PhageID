# DefenseFinder
rule defenseFinder_lastversion:
    output: os.path.join(RESULTS_DIR, "logs", "defensefinder_lastversion.txt")
    log: os.path.join(RESULTS_DIR, "logs", "defensefinder_lastversion.log")
    conda: os.path.join(ENV_DIR, "defensefinder.yaml")
    shell:
        """(date && defense-finder update && date) > {output}"""

rule antiDefenseFinder:
    output: os.path.join(RESULTS_DIR, "{sample}", "defenseFinder", "phanotate_defense_finder_systems.tsv")
    input: 
        pharokka = rules.pharokka_phage.output.faa,
        defenseFinder_lastversion = rules.defenseFinder_lastversion.output
    log: os.path.join(RESULTS_DIR, "logs", "{sample}_defenseFinder.log")
    conda: os.path.join(ENV_DIR, "defensefinder.yaml")
    threads: 4    
    shell:
        """(date && 
        defense-finder run -w {threads} --preserve-raw --antidefensefinder -o $(dirname {output}) {input.pharokka} &&
        date) &> {log}"""
