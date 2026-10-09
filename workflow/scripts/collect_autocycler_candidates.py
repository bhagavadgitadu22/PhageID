#!/usr/bin/env python3
"""Collect successful per-job Autocycler candidates and record all outcomes."""
import argparse
import csv
import shutil
from pathlib import Path


def collect(candidates, output, allow_insufficient=False):
    output = Path(output)
    assemblies = output / 'assemblies'
    assemblies.mkdir(parents=True, exist_ok=True)
    # Clear only generated collection files when rebuilding after a candidate retry.
    for previous in assemblies.glob('*.fasta'):
        previous.unlink()
    rows = []
    for directory in candidates:
        directory = Path(directory)
        status = (directory / 'status.txt').read_text().strip()
        if status not in {'success', 'failed'}:
            raise ValueError(f'Invalid candidate status: {directory}')
        name = directory.parent.name + '_' + directory.name
        if status == 'success':
            fasta = directory / 'assembly.fasta'
            if not fasta.is_file() or fasta.stat().st_size == 0:
                raise ValueError(f'Successful candidate has no assembly: {directory}')
            shutil.copyfile(fasta, assemblies / (name + '.fasta'))
        rows.append({'Assembler': directory.parent.name, 'Subset': directory.name,
                     'Status': status, 'Candidate directory': str(directory)})
    with open(output / 'candidate_status.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['Assembler', 'Subset', 'Status', 'Candidate directory'],
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    successes = sum(row['Status'] == 'success' for row in rows)
    if successes < 2 and not allow_insufficient:
        raise ValueError(f'Autocycler requires multiple successful input assemblies; got {successes}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--allow-insufficient', action='store_true')
    args = parser.parse_args()
    collect(args.candidates, args.output, args.allow_insufficient)


if __name__ == '__main__':
    main()
