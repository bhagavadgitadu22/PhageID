"""Run theBIGbam enrichment, passing through annotations when a CSV has no hits."""

import argparse
import csv
import re
import shutil
import subprocess
from pathlib import Path


def restore_annotation_pipes(path, attribute):
    """Keep pipe-separated predictions readable in one GFF3 attribute."""
    gff = Path(path)
    lines = gff.read_text().splitlines(keepends=True)
    in_fasta = False
    for index, line in enumerate(lines):
        if line.startswith("##FASTA"):
            in_fasta = True
        if in_fasta or line.startswith("#"):
            continue
        fields = line.split("\t", 8)
        if len(fields) != 9:
            continue
        attrs = fields[8].split(";")
        for attr_index, item in enumerate(attrs):
            if item.startswith(attribute + "="):
                attrs[attr_index] = re.sub(r"%7c", "|", item, flags=re.IGNORECASE)
        fields[8] = ";".join(attrs)
        lines[index] = "\t".join(fields)
    gff.write_text("".join(lines))


def enrich(annotations, csv_path, output, prefix, literal_pipes=False):
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        next(reader, None)  # header
        has_rows = any(any(value.strip() for value in row) for row in reader)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    if not has_rows:
        print(f"No annotation rows in {csv_path}; copying input annotations unchanged.", flush=True)
        shutil.copyfile(annotations, output)
        return
    subprocess.run([
        "thebigbam", "add-contig-annotations", "-g", str(annotations),
        "--csv", str(csv_path), "--match-by", "feature_type,locus_tag",
        "--prefix", prefix, "--keep-multiple", "-o", str(output),
    ], check=True)
    if literal_pipes:
        restore_annotation_pipes(output, prefix + "Annotation")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--literal-pipes", action="store_true",
                        help="Restore literal | in the prefixed Annotation GFF attribute")
    args = parser.parse_args()
    enrich(args.annotations, args.csv, args.output, args.prefix, args.literal_pipes)
