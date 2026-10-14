#!/usr/bin/env python3
"""Prepare one sequence-containing GFF per Pharokka contig for LoVis4u."""

import argparse
import hashlib
import re
import shutil
import tempfile
from pathlib import Path


def safe_label(value):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)[:48]


def split_gff(source, sample, destination):
    features = {}
    sequences = {}
    current_id = None
    in_fasta = False

    with source.open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                if in_fasta:
                    raise ValueError(f"{source}:{number}: duplicate ##FASTA section")
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith(">"):
                    current_id = line[1:].split(maxsplit=1)[0]
                    if not current_id or current_id in sequences:
                        raise ValueError(f"{source}:{number}: empty or duplicate FASTA ID {current_id!r}")
                    sequences[current_id] = [line]
                elif current_id is None:
                    raise ValueError(f"{source}:{number}: FASTA sequence without a header")
                else:
                    sequences[current_id].append(line)
            elif line.strip() and not line.startswith("#"):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 9:
                    raise ValueError(f"{source}:{number}: expected at least nine GFF fields")
                features.setdefault(fields[0], []).append((number, line, fields))

    if not sequences:
        raise ValueError(f"{source}: missing embedded FASTA sequences")
    missing = set(features) - set(sequences)
    if missing:
        raise ValueError(f"{source}: missing FASTA sequence for {', '.join(sorted(missing))}")

    for contig, fasta_lines in sequences.items():
        length = sum(len(line.strip()) for line in fasta_lines[1:])
        if not length:
            raise ValueError(f"{source}: empty FASTA sequence for {contig}")
        for number, _, fields in features.get(contig, []):
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{source}:{number}: invalid feature coordinates") from error
            if not 1 <= start <= end <= length:
                raise ValueError(f"{source}:{number}: feature exceeds {contig} length {length}")

        if not features.get(contig):
            continue  # LoVis4u skips loci without annotation features.

        digest = hashlib.sha256(f"{sample}\0{contig}".encode()).hexdigest()[:12]
        name = f"{safe_label(sample)}__{safe_label(contig)}__{digest}.gff"
        with (destination / name).open("w") as output:
            output.write("##gff-version 3\n")
            output.write(f"##sequence-region {contig} 1 {length}\n")
            for _, line, _ in features.get(contig, []):
                output.write(line)
                if not line.endswith("\n"):
                    output.write("\n")
            output.write("##FASTA\n")
            output.writelines(fasta_lines)
            if not fasta_lines[-1].endswith("\n"):
                output.write("\n")
    return sum(bool(features.get(contig)) for contig in sequences)


def prepare(gff_paths, output_dir):
    if not gff_paths:
        raise ValueError("No Pharokka GFF inputs supplied")
    output_dir = Path(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    samples = set()
    with tempfile.TemporaryDirectory(prefix="lovis4u_gff_", dir=output_dir.parent) as temporary:
        staging = Path(temporary)
        count = 0
        for source in map(Path, gff_paths):
            sample = source.parent.parent.name
            if sample in samples:
                raise ValueError(f"Duplicate Pharokka sample: {sample}")
            samples.add(sample)
            count += split_gff(source, sample, staging)
        if not count:
            raise ValueError("No annotated contigs found in Pharokka GFF inputs")
        if output_dir.is_symlink():
            raise ValueError(f"Refusing to replace symlink: {output_dir}")
        if output_dir.exists():
            shutil.rmtree(output_dir)
        staging.rename(output_dir)
    print(f"Prepared {count} single-contig GFF files from {len(samples)} samples")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("gff", nargs="+")
    args = parser.parse_args()
    prepare(args.gff, args.output_dir)


if __name__ == "__main__":
    main()
