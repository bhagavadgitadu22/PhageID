#!/usr/bin/env python3
"""Prepare representative FASTAs and sample-to-representative mapping jobs."""
import argparse
import csv
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO


def prepare(representatives_fasta, clusters, sample_fastas, output_dir):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(representatives_fasta, 'fasta'))
    representatives = {record.id: record for record in records}
    if len(representatives) != len(records):
        raise ValueError('Duplicate representative IDs')
    owners = {}
    for sample, fasta in sample_fastas:
        for record in SeqIO.parse(fasta, 'fasta'):
            if record.id in owners:
                raise ValueError(f'Duplicate input contig: {record.id}')
            owners[record.id] = sample
    assignments = {}
    pairs = set()
    seen_representatives = set()
    with open(clusters) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2:
                raise ValueError(f'Invalid cluster row {line_number}')
            representative, members = fields
            if representative not in representatives or representative in seen_representatives:
                raise ValueError(f'Unknown or duplicate representative: {representative}')
            seen_representatives.add(representative)
            member_ids = members.split(',')
            if representative not in member_ids:
                raise ValueError(f'Representative missing from its cluster: {representative}')
            for member in member_ids:
                if member not in owners or member in assignments:
                    raise ValueError(f'Unknown or duplicated cluster member: {member}')
                assignments[member] = representative
                pairs.add((owners[member], representative))
    if set(assignments) != set(owners) or seen_representatives != set(representatives):
        raise ValueError('Cluster membership does not cover the input contigs and representatives')
    reference_dir = output / 'references'
    reference_dir.mkdir(exist_ok=True)
    for representative, record in representatives.items():
        SeqIO.write([record], reference_dir / (quote(representative, safe='_-.') + '.fna'), 'fasta')
    with (output / 'sample_representatives.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'representative', 'representative_id'])
        for sample, representative in sorted(pairs):
            writer.writerow([sample, quote(representative, safe='_-.'), representative])


def read_pairs(manifest):
    with open(manifest, newline='') as handle:
        return [(row['sample'], row['representative']) for row in csv.DictReader(handle, delimiter='\t')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--representatives-fasta', required=True)
    parser.add_argument('--clusters', required=True)
    parser.add_argument('--sample-fasta', action='append', nargs=2, metavar=('SAMPLE', 'FASTA'), default=[])
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    prepare(args.representatives_fasta, args.clusters, args.sample_fasta, args.output_dir)


if __name__ == '__main__':
    main()
