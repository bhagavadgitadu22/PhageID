#!/usr/bin/env python3
"""Report Pharokka function counts and CDS coding density by viral contig."""

import argparse
import csv
from pathlib import Path


CATEGORIES = [
    "CDS",
    "connector",
    "DNA, RNA and nucleotide metabolism",
    "head and packaging",
    "integration and excision",
    "lysis",
    "moron, auxiliary metabolic gene and host takeover",
    "other",
    "tail",
    "transcription regulation",
    "unknown function",
    "tRNAs",
    "CRISPRs",
    "tmRNAs",
    "VFDB_Virulence_Factors",
    "CARD_AMR_Genes",
]
FIELDS = ["Viral_contig", "Sample", "Coding_density_percent", *CATEGORIES]


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_lengths(path):
    lengths = {}
    contig = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                contig = line[1:].split()[0]
                if not contig or contig in lengths:
                    raise ValueError(f"Missing or duplicate contig ID in {path}: {contig}")
                lengths[contig] = 0
            elif line:
                if contig is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                lengths[contig] += len(line)
    if any(length == 0 for length in lengths.values()):
        raise ValueError(f"Zero-length contig in {path}")
    return lengths


def cds_intervals(path, lengths):
    intervals = {contig: [] for contig in lengths}
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid GFF row at {path}:{number}")
            if fields[2] != "CDS":
                continue
            if len(fields) < 9 or any(fields[9:]):
                raise ValueError(f"Invalid CDS GFF row at {path}:{number}")
            contig = fields[0]
            if contig not in lengths:
                raise ValueError(f"CDS belongs to unknown contig {contig} at {path}:{number}")
            start, end = int(fields[3]), int(fields[4])
            if not 1 <= start <= end <= lengths[contig]:
                raise ValueError(f"CDS coordinates outside {contig} at {path}:{number}")
            intervals[contig].append((start, end))
    return intervals


def covered_bases(intervals):
    covered = 0
    last_end = 0
    for start, end in sorted(intervals):
        covered += max(0, end - max(start, last_end + 1) + 1)
        last_end = max(last_end, end)
    return covered


def function_counts(path, lengths):
    counts = {contig: dict.fromkeys(CATEGORIES, 0) for contig in lengths}
    seen = set()
    for row in read_tsv(path):
        category, contig = row["Description"], row["contig"]
        if contig not in counts:
            raise ValueError(f"Function count belongs to unknown contig {contig} in {path}")
        if category not in CATEGORIES:
            raise ValueError(f"Unknown Pharokka function category {category} in {path}")
        if (contig, category) in seen:
            raise ValueError(f"Duplicate Pharokka category {category} for {contig} in {path}")
        value = int(row["Count"])
        if value < 0:
            raise ValueError(f"Negative Pharokka count for {contig} in {path}")
        counts[contig][category] = value
        seen.add((contig, category))
    return counts


def sample_report(sample, fasta, gff, functions, output):
    lengths = fasta_lengths(fasta)
    intervals = cds_intervals(gff, lengths)
    counts = function_counts(functions, lengths)
    rows = []
    for contig, length in lengths.items():
        density = 100 * covered_bases(intervals[contig]) / length
        rows.append({"Viral_contig": contig, "Sample": sample,
                     "Coding_density_percent": f"{density:.2f}", **counts[contig]})
    write_tsv(output, rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sample", "fasta", "gff", "functions", "output"):
        parser.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    sample_report(args.sample, args.fasta, args.gff, args.functions, args.output)


if __name__ == "__main__":
    main()
