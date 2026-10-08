#!/usr/bin/env python3
"""Left-join Empathi and Sublyme CSV predictions using their first two columns."""

import argparse
import csv
import logging
from collections import Counter
from pathlib import Path

LYSIS = {"holin", "lysin", "endolysin", "spanin", "val", "lysis_inhibitor"}
DEFENSE = {"anti-restriction", "toxin", "super_infection", "crispr", "sir2"}

def tokens(annotation):
    # Empty CSV fields and common serialized missing values are not labels.
    return [label for part in annotation.split("|")
            if (label := part.strip()) and label.lower() not in {"na", "nan", "none", "null"}]

def merge_annotations(empathi, sublyme):
    labels = list(dict.fromkeys(
        [label for label in tokens(empathi) if label != "endolysin"] + tokens(sublyme)
    ))
    if any(label != "unknown" for label in labels):
        labels = [label for label in labels if label != "unknown"]
    for members, category in ((LYSIS, "lysis"), (DEFENSE, "defense_systems")):
        if members.intersection(labels) and category not in labels:
            labels.append(category)
    return "|".join(dict.fromkeys(label.replace("_", " ") for label in labels))

def read_predictions(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        if len(header) < 2:
            raise ValueError(f"{path}: expected a header with at least two columns")
        for row in reader:
            if not row:
                continue
            if len(row) != len(header) or not row[0].strip():
                raise ValueError(f"{path}: malformed prediction at line {reader.line_num}")
            yield row[0].strip(), row[1]

def combine(empathi_path, sublyme_path, output_path):
    sublyme = {}
    counts = Counter()
    for protein, annotation in read_predictions(sublyme_path):
        counts[protein] += 1
        sublyme.setdefault(protein, []).extend(tokens(annotation))
    for protein, count in counts.items():
        if count > 1:
            logging.warning("Duplicate Sublyme protein ID %r: %d rows; combining labels in file order", protein, count)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    written = unmatched = 0
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Protein", "Annotation"])
        for protein, annotation in read_predictions(empathi_path):
            unmatched += protein not in sublyme
            writer.writerow([protein, merge_annotations(annotation, "|".join(sublyme.get(protein, [])))])
            written += 1
    logging.info("Wrote %d Empathi rows; %d had no Sublyme match", written, unmatched)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--empathi", required=True)
    parser.add_argument("--sublyme", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    combine(args.empathi, args.sublyme, args.output)

if __name__ == "__main__":
    main()
