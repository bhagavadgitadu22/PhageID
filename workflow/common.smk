import os


def sample_path(value):
    import os

    value = value.strip()
    if value.lower() in {"", "na", "none", "-"}:
        return ""
    if os.path.isabs(value):
        return value
    root = globals().get("SAMPLE_DATA_ROOT") or os.path.dirname(
        os.path.abspath(SAMPLES_FILE)
    )
    return os.path.normpath(os.path.join(root, value))


def read_input_files(path):
    if not path:
        raise ValueError("No read input supplied")
    if not os.path.isdir(path):
        return [path]
    files = [
        os.path.join(path, name)
        for name in sorted(os.listdir(path))
        if name.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))
        and os.path.isfile(os.path.join(path, name))
    ]
    if not files:
        raise ValueError(f"No FASTQ files found in {path}")
    return files


def sample_reads(wildcards):
    pattern = (
        rules.preprocess_reads_porechop.output[0]
        if READ_FILES[wildcards.sample]
        else rules.preprocess_reads_fastp.output.reads
    )
    return pattern.format(sample=wildcards.sample)


def assembly_reads(wc):
    if wc.readset == "filtered":
        return rules.remove_bacterial_contamination.output[0].format(sample=wc.sample)
    return sample_reads(wc)


def candidate_assembly(wc):
    if wc.assembler in {"flye", "autocycler"} and not READ_FILES[wc.sample]:
        raise ValueError(f"{wc.assembler} requires long reads for {wc.sample}")
    return os.path.join(
        RESULTS_DIR,
        wc.sample,
        "assembly_attempts",
        wc.readset,
        wc.assembler,
        "assembly.fasta",
    )


def candidate_choice(sample, readset, dependencies=None):
    import csv

    primary = "flye" if READ_FILES[sample] else "spades"

    def quality(assembler):
        candidate = checkpoints.candidate_viral_quality.get(
            sample=sample, readset=readset, assembler=assembler
        ).output
        if dependencies is not None:
            dependencies.extend([str(candidate.quality), str(candidate.status)])
        with open(candidate.status) as handle:
            status = handle.read().strip()
        if status == "failed":
            return -1
        if status != "success":
            raise ValueError(f"Invalid assembly status: {status}")
        with open(candidate.quality) as handle:
            return max(
                (
                    CHECKV_QUALITY_RANK.get(row["checkv_quality"], 1)
                    for row in csv.DictReader(handle, delimiter="\t")
                ),
                default=0,
            )

    primary_quality = quality(primary)
    if primary == "spades" or primary_quality >= CHECKV_QUALITY_RANK["High-quality"]:
        return primary, primary_quality
    autocycler_quality = quality("autocycler")
    # Prefer the earlier attempt when the best categories are equal.
    if autocycler_quality > primary_quality:
        return "autocycler", autocycler_quality
    return primary, primary_quality


def selected_candidate_inputs(wc, dependencies=None):
    import csv

    filter_report = checkpoints.assembly_read_status.get(sample=wc.sample).output[0]
    if dependencies is not None:
        dependencies.append(str(filter_report))
    with open(filter_report) as handle:
        stats = next(csv.DictReader(handle, delimiter="\t"))
    readset = (
        "filtered"
        if HOSTS_LIST[wc.sample] and int(stats["Bacterial free reads"]) > 0
        else "all"
    )
    assembler, best_quality = candidate_choice(wc.sample, readset, dependencies)
    if readset == "filtered" and best_quality < CHECKV_QUALITY_RANK["High-quality"]:
        all_assembler, all_quality = candidate_choice(wc.sample, "all", dependencies)
        if all_quality > best_quality:
            readset, assembler = "all", all_assembler
    values = dict(rules.fix_circular_viral_contigs_per_sample_candidate.output.items())
    values.update(dict(rules.genomad_candidate.output.items()))
    values.update(dict(rules.checkv_candidate.output.items()))
    values["status"] = os.path.join(
        RESULTS_DIR,
        "{sample}",
        "assembly_attempts",
        "{readset}",
        "{assembler}",
        "assembly_status.txt",
    )
    values["assembly"] = os.path.join(
        RESULTS_DIR,
        "{sample}",
        "assembly_attempts",
        "{readset}",
        "{assembler}",
        "assembly.fasta",
    )
    if assembler == "flye":
        values["flye_info"] = rules.assembly_reads_flye.output.info
        values["used_reads"] = rules.subsample_reads_flye.output.reads
    else:
        values["used_reads"] = (
            rules.remove_bacterial_contamination.output[0]
            if readset == "filtered"
            else sample_reads(wc)
        )
    values["filter_report"] = filter_report
    return {
        name: str(path).format(sample=wc.sample, assembler=assembler, readset=readset)
        for name, path in values.items()
    }


def selection_manifest_inputs(wc):
    dependencies = []
    sources = selected_candidate_inputs(wc, dependencies)
    # Track every assessed candidate, including those that lost the comparison.
    return list(dict.fromkeys(dependencies + list(sources.values())))


def sample_mapper(wildcards):
    return "minimap2-ont" if READ_FILES[wildcards.sample] else "minimap2-sr-secondary"


def manual_viral_fasta(wildcards):
    path = os.path.join(
        RESULTS_DIR, wildcards.sample, "manual_corrections", "circular_viruses.fasta"
    )
    return [path] if os.path.isfile(path) else []


def active_viral_samples(wildcards):
    import csv

    sample_report = checkpoints.viral_report_global.get().output.samples
    with open(sample_report, newline="") as handle:
        return [
            row["Sample"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["Status"] == "ready"
        ]


def active_sample_files(pattern):
    return lambda wildcards: expand(pattern, sample=active_viral_samples(wildcards))


def per_sample_pipeline_inputs(wildcards):
    sample_report = checkpoints.viral_report_global.get().output.samples
    if wildcards.sample not in PHAGES_LIST:
        raise ValueError(f"Unknown sample: {wildcards.sample}")
    if wildcards.sample not in active_viral_samples(wildcards):
        return [sample_report]
    outputs = [
        *rules.phageterm.output,
        *rules.viral_annotations_report.output,
        *rules.pharokka_plot.output,
        *rules.genotate_pharokka_plot.output,
        *rules.combined_viral_annotations_per_sample.output,
    ]
    return expand(outputs, sample=wildcards.sample)


def comparison_pipeline_inputs(wildcards):
    sample_report = checkpoints.viral_report_global.get().output.samples
    if not active_viral_samples(wildcards):
        return [sample_report]
    return [
        sample_report,
        rules.ani_for_dereplication.output.ani_results,
        *rules.viruses_dereplicated.output,
        *rules.dereplication_report.output,
        *rules.thebigbam_calculate.output,
        *rules.thebigbam_calculate_post_dereplication.output,
        *rules.lovis4u_pharokka.output,
        *rules.lovis4u_pharokka_post_dereplication.output,
        *rules.vcontact3.output,
    ]


def assembly_databases():
    # Plassembler is only used by the long-read Autocycler rescue route.
    return ASSEMBLY_DATABASES + (
        ("plassembler_db",) if any(READ_FILES.values()) else ()
    )


def validate_inputs(required_databases):
    problems = []
    for sample in PHAGES_LIST:
        for label, path in (
            ("long reads", READ_FILES[sample]),
            ("short reads", SHORT_READS[sample]),
        ):
            if not path:
                continue
            if os.path.isdir(path):
                try:
                    read_input_files(path)
                except ValueError as error:
                    problems.append(f"{sample} {label}: {error}")
            elif not os.path.isfile(path):
                problems.append(f"{sample} {label}: missing file or directory: {path}")
        host = HOSTS_LIST[sample]
        if host and not os.path.isfile(host):
            problems.append(f"{sample} host reference: missing file: {host}")

    for name in required_databases:
        path = config[name]
        if not os.path.isdir(path):
            problems.append(f"{name}: missing database directory: {path}")
    if "genomad_db" in required_databases:
        genomad_marker = os.path.join(
            config["genomad_db"], "genomad_marker_metadata.tsv"
        )
        if os.path.isdir(config["genomad_db"]) and not os.path.isfile(genomad_marker):
            problems.append(f"genomad_db: missing marker metadata: {genomad_marker}")
    if problems:
        raise ValueError(
            "Missing pipeline inputs:\n" + "\n".join(f"- {item}" for item in problems)
        )


def preflight_inputs(required_databases):
    def check(_wildcards):
        validate_inputs(required_databases)
        return []

    return check


def autocycler_candidate_jobs(wc):
    subsets = int(rules.autocycler_subsample.params.subset_count)
    if subsets < 1:
        raise ValueError("subset_count must be at least 1")
    return expand(
        rules.autocycler_candidate_assembly.output[0],
        sample=wc.sample,
        readset=wc.readset,
        candidate_assembler=[
            "canu",
            "flye",
            "metamdbg",
            "miniasm",
            "necat",
            "nextdenovo",
            "plassembler",
            "raven",
        ],
        subset=[f"{i:02d}" for i in range(1, subsets + 1)],
    )


def post_dereplication_pairs(wildcards):
    import csv

    prepared = checkpoints.prepare_post_dereplication_mapping.get().output.prepared
    with open(
        os.path.join(prepared, "sample_representatives.tsv"), newline=""
    ) as handle:
        return [
            (row["sample"], row["representative"])
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def post_dereplication_bams(wildcards):
    return [
        os.path.join(
            RESULTS_DIR,
            "minimap2",
            "thebigbam_post_dereplication",
            f"{sample}_on_{representative}.bam",
        )
        for sample, representative in post_dereplication_pairs(wildcards)
    ]
