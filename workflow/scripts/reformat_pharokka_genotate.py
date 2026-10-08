#!/usr/bin/env python3
"""Convert Genotate/Pharokka protein annotations to a multi-contig GFF3."""
import argparse
import csv
import re
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO

PATTERN = re.compile(r"(.*?)_CDS_\[(complement\()?(\d+)\.\.(\d+)(\))?\]$")


def convert(genes, fasta, output):
    records = list(SeqIO.parse(fasta, "fasta"))
    lengths = {record.id: len(record.seq) for record in records}
    if len(lengths) != len(records):
        raise ValueError("Duplicate FASTA contig IDs")
    features = {name: [] for name in lengths}
    with open(genes, newline="") as handle:
        for line, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            match = PATTERN.fullmatch(row[0])
            if not match:
                if line == 1 and row[0].lower() in {"name", "id", "protein", "protein_id", "gene", "gene_id"}:
                    continue
                raise ValueError(f"{genes}:{line}: unrecognized protein ID {row[0]!r}")
            if len(row) < 5:
                raise ValueError(f"{genes}:{line}: expected at least five columns")
            contig, complement, start, end, closing = match.groups()
            if bool(complement) != bool(closing):
                raise ValueError(f"{genes}:{line}: malformed complement coordinates")
            start, end = int(start), int(end)
            if contig not in lengths or not 1 <= start <= end <= lengths[contig]:
                raise ValueError(f"{genes}:{line}: unknown contig or invalid coordinates")
            category = row[4].replace("DNA, RNA and nucleotide metabolism", "DNA").replace("moron, auxiliary metabolic gene and host takeover", "moron")
            attributes = {"ID": row[0], "phrog": row[2], "function": category, "product": row[3]}
            attrs = ";".join(f"{key}={quote(value, safe=':_-.')}" for key, value in attributes.items())
            # Each Genotate prediction is a complete CDS; genomic reading frame
            # belongs in plotting logic, not the GFF phase (bases to skip).
            features[contig].append("\t".join([contig, "genotate_0.15", "CDS", str(start), str(end), ".", "-" if complement else "+", "0", attrs]))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        handle.write("##gff-version 3\n")
        for record in records:
            handle.write(f"##sequence-region {record.id} 1 {len(record.seq)}\n")
        for record in records:
            for feature in features[record.id]:
                handle.write(feature + "\n")
        handle.write("##FASTA\n")
        SeqIO.write(records, handle, "fasta")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    convert(args.genes, args.fasta, args.output)


if __name__ == "__main__":
    main()
