#!/usr/bin/env python3
"""Publish an optional manual viral FASTA and record every sequence change."""

import argparse
import csv
import shutil
from pathlib import Path


REPORT_FIELDS = ["Viral_contig", "Auto_corrected_bp", "Final_bp",
                 "Manual_correction", "Manual_length_change_bp"]


def fasta_sequences(path):
    sequences = {}
    name = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name or name in sequences:
                    raise ValueError(f"Missing or duplicate FASTA ID in {path}: {name}")
                sequences[name] = ""
            elif name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                sequences[name] += "".join(line.split())
    for contig, sequence in sequences.items():
        if not sequence:
            raise ValueError(f"Empty sequence for contig {contig} in {path}")
    return sequences


def correction_rows(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    names = [row["contig_id"] for row in rows]
    if len(set(names)) != len(names):
        raise ValueError(f"Duplicate contig ID in {path}")
    return {row["contig_id"]: row for row in rows}


def publish(automatic, concatemer, dtr, manual, output, concatemer_output,
            dtr_output, report):
    automatic_sequences = fasta_sequences(automatic)
    final_sequences = fasta_sequences(manual) if manual else automatic_sequences
    if list(final_sequences) != list(automatic_sequences):
        raise ValueError("Manual FASTA must have exactly the automatic contig IDs in the same order")
    concatemers = correction_rows(concatemer)
    dtrs = correction_rows(dtr)
    if set(concatemers) != set(automatic_sequences) or set(dtrs) != set(automatic_sequences):
        raise ValueError("Correction report IDs do not match automatic FASTA")
    rows = []
    for name, sequence in automatic_sequences.items():
        if int(concatemers[name]["corrected_length"]) != int(dtrs[name]["original_length"]):
            raise ValueError(f"Correction stages disagree for {name}")
        if int(dtrs[name]["corrected_length"]) != len(sequence):
            raise ValueError(f"Automatic corrected length differs from DTR report for {name}")
        final = final_sequences[name]
        rows.append([name, len(sequence), len(final),
                     "Yes" if final != sequence else "No", len(final) - len(sequence)])
    for destination in (output, concatemer_output, dtr_output, report):
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(manual or automatic, output)
    shutil.copyfile(concatemer, concatemer_output)
    shutil.copyfile(dtr, dtr_output)
    with open(report, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(REPORT_FIELDS)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("automatic", "concatemer", "dtr", "output", "concatemer-output",
                     "dtr-output", "report"):
        parser.add_argument("--" + argument, required=True)
    parser.add_argument("--manual")
    args = parser.parse_args()
    publish(args.automatic, args.concatemer, args.dtr, args.manual, args.output,
            args.concatemer_output, args.dtr_output, args.report)


if __name__ == "__main__":
    main()
