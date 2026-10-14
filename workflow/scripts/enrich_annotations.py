"""Run theBIGbam enrichment, passing through annotations when a CSV has no hits."""

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def enrich(annotations, csv_path, output, prefix):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    enrich(args.annotations, args.csv, args.output, args.prefix)
