#!/usr/bin/env python3
"""Convert a delimited protein annotation table for theBIGbam enrichment."""

import argparse
import csv
import logging
import re
from pathlib import Path

ID_COLUMNS = (
    "locus_tag",
    "hit_id",
    "protein_id",
    "protein",
    "gene_id",
    "gene",
    "query",
    "sequence_id",
    "seq_id",
    "name",
    "id",
    "",
    "Unnamed: 0",
)

def delimiter_for(path):
    with open(path, newline="") as handle:
        sample = handle.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        return "\t" if Path(path).suffix.lower() in {".tsv", ".tab"} else ","

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--id-column", help="Column containing the Pharokka/Bakta protein ID (auto-detected by default)")
    parser.add_argument("--checkamg-pharokka-fasta", help="Map raw CheckAMG IDs to actual Pharokka FASTA IDs and retain checkamg_id")
    args = parser.parse_args()

    pharokka_ids = {}
    if args.checkamg_pharokka_fasta:
        with open(args.checkamg_pharokka_fasta) as fasta:
            for line in fasta:
                if not line.startswith(">"):
                    continue
                protein_id = line[1:].split()[0]
                match = re.fullmatch(r"(.+)_CDS_(\d+)", protein_id)
                if not match:
                    raise ValueError(f"Unexpected Pharokka protein ID: {protein_id!r}")
                raw_id = f"{match[1]}_{int(match[2])}"
                if raw_id in pharokka_ids:
                    raise ValueError(f"Ambiguous Pharokka mapping for {raw_id!r}")
                pharokka_ids[raw_id] = protein_id
        if not pharokka_ids:
            raise ValueError("No Pharokka protein IDs found in FASTA")

    logging.info("Reading feature annotations from %s", args.input)
    with open(args.input, newline="") as source:
        reader = csv.DictReader(source, delimiter=delimiter_for(args.input))
        columns = list(reader.fieldnames or [])
        if not columns:
            raise ValueError(f"No header found in {args.input}")

        if args.id_column:
            id_column = args.id_column
            if id_column not in columns:
                raise ValueError(
                    f"ID column {id_column!r} is absent from {args.input}; found {columns}"
                )
        else:
            columns_by_lowercase = {name.lower(): name for name in columns}
            id_column = next(
                (
                    columns_by_lowercase[name.lower()]
                    for name in ID_COLUMNS
                    if name.lower() in columns_by_lowercase
                ),
                None,
            )
            if id_column is None:
                raise ValueError(
                    f"Could not identify a protein ID column in {args.input}; found {columns}"
                )

        # Rename the ID column and replace any existing feature metadata.
        annotation_columns = [
            column for column in columns
            if column not in {id_column, "feature_type", "locus_tag"}
        ]
        if args.checkamg_pharokka_fasta and "checkamg_id" in annotation_columns:
            raise ValueError("Input already contains a checkamg_id annotation column")
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        written = skipped = 0
        with output.open("w", newline="") as target:
            writer = csv.writer(target, lineterminator="\n")
            original_id_columns = ["checkamg_id"] if args.checkamg_pharokka_fasta else []
            writer.writerow(["feature_type", "locus_tag", *original_id_columns, *annotation_columns])
            for row in reader:
                locus_tag = (row.get(id_column) or "").strip()
                if not locus_tag:
                    skipped += 1
                    continue
                original_id_values = []
                if args.checkamg_pharokka_fasta:
                    original_id_values = [locus_tag]
                    if locus_tag not in pharokka_ids:
                        raise ValueError(
                            f"{args.input}:{reader.line_num}: no Pharokka protein for CheckAMG ID {locus_tag!r}"
                        )
                    locus_tag = pharokka_ids[locus_tag]
                writer.writerow([
                    "CDS", locus_tag, *original_id_values,
                    *(row.get(column, "") for column in annotation_columns),
                ])
                written += 1
        logging.info(
            "Wrote %s features to %s using ID column %r; skipped %s rows without IDs",
            written, output, id_column, skipped,
        )

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        main()
    except Exception:
        logging.exception("Feature annotation conversion failed")
        raise SystemExit(1)
