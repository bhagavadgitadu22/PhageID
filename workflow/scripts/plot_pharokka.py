#!/usr/bin/env python3
"""Plot Pharokka contigs that have a labelable CDS product."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

from Bio import SeqIO


def has_labelable_cds(record):
    return any(
        feature.type == "CDS"
        and any(
            product.strip() and not product.strip().lower().startswith(("hypothetical", "unknown"))
            for product in feature.qualifiers.get("product", [])
        )
        for feature in record.features
    )


def plot(genbank, output_dir):
    records = list(SeqIO.parse(genbank, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records in {genbank}")
    selected = [record for record in records if has_labelable_cds(record)]
    output = Path(output_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pharokka_plot_", dir=output.parent) as staging:
        plots = Path(staging) / "plots"
        if selected:
            filtered = Path(staging) / "labelable.gbk"
            SeqIO.write(selected, filtered, "genbank")
            subprocess.run(["pharokka", "multiplot", "-g", str(filtered), "-o", str(plots)], check=True)
            if not plots.is_dir():
                raise RuntimeError("Pharokka succeeded without creating a plots directory")
        else:
            plots.mkdir()
        if output.exists():
            shutil.rmtree(output)
        shutil.move(str(plots), output)
    print(f"Plotted {len(selected)} of {len(records)} contigs; skipped {len(records) - len(selected)} without named CDS products")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genbank", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    plot(args.genbank, args.output_dir)


if __name__ == "__main__":
    main()
