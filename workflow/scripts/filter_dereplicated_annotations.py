#!/usr/bin/env python3
"""Merge enriched GFF features for contigs selected by one or more FASTAs."""
import argparse
import os
import tempfile
from pathlib import Path

from Bio import SeqIO


def filter_annotations(annotations, representatives_fasta, output):
    fasta_paths = [representatives_fasta] if isinstance(representatives_fasta, (str, Path)) else representatives_fasta
    records = [record for path in fasta_paths for record in SeqIO.parse(path, 'fasta')]
    lengths = {record.id: len(record.seq) for record in records}
    if len(lengths) != len(records):
        raise ValueError('Duplicate selected FASTA IDs')
    features = {record.id: [] for record in records}
    sources = {}
    for annotation in annotations:
        with open(annotation) as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith('##FASTA'):
                    break  # Emit sequences from the authoritative selected FASTAs.
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 9:
                    raise ValueError(f'{annotation}:{line_number}: expected nine GFF columns')
                contig = fields[0]
                if contig not in lengths:
                    continue
                start, end = int(fields[3]), int(fields[4])
                if not 1 <= start <= end <= lengths[contig]:
                    raise ValueError(f'{annotation}:{line_number}: coordinates outside representative {contig}')
                if contig in sources and sources[contig] != str(annotation):
                    raise ValueError(f'Representative {contig} occurs in multiple annotation files')
                sources[contig] = str(annotation)
                features[contig].append('\t'.join(fields) + '\n')
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write('##gff-version 3\n')
            for record in records:
                handle.write(f'##sequence-region {record.id} 1 {len(record.seq)}\n')
            for record in records:
                handle.writelines(features[record.id])
            if records:
                handle.write('##FASTA\n')
                SeqIO.write(records, handle, 'fasta')
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', nargs='*', default=[])
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--representatives-fasta')
    selection.add_argument('--fasta', nargs='+')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    filter_annotations(args.annotations, args.fasta or args.representatives_fasta, args.output)


if __name__ == '__main__':
    main()
