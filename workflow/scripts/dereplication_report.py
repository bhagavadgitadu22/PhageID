#!/usr/bin/env python3
"""Report cluster membership and representative selection for every viral contig."""
import argparse
import csv
from pathlib import Path
from viral_contig_report import fasta_lengths

FIELDS = ['Viral contig', 'Length', 'Selected', 'Representative',
          'Representative length', 'Cluster size', 'ANI identity (%)',
          'Aligned fraction (%)', 'Representative aligned fraction (%)']


def write_report(fasta, clusters, representatives, ani, output):
    lengths = fasta_lengths(fasta)
    selected = fasta_lengths(representatives)
    membership = {}
    cluster_sizes = {}
    with open(clusters) as handle:
        for line in handle:
            if not line.strip():
                continue
            representative, members = line.rstrip('\r\n').split('\t')
            members = members.split(',')
            if representative in cluster_sizes or representative not in members:
                raise ValueError(f'Invalid cluster representative: {representative}')
            cluster_sizes[representative] = len(members)
            for member in members:
                if member not in lengths or member in membership:
                    raise ValueError(f'Unknown or duplicated cluster member: {member}')
                membership[member] = representative
    if set(membership) != set(lengths):
        raise ValueError('Clusters do not cover all input contigs')
    if set(selected) != set(cluster_sizes):
        raise ValueError('Selected contigs do not match cluster representatives')
    if any(length != lengths[name] for name, length in selected.items()):
        raise ValueError('Representative lengths differ from the input FASTA')
    # Clustering uses directed representative (query) -> member (target) edges.
    alignments = {}
    with open(ani, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            query, target = row['qname'], row['tname']
            if membership.get(target) != query or query == target:
                continue
            key = (query, target)
            if key in alignments:
                raise ValueError(f'Duplicate ANI pair: {query}, {target}')
            alignments[key] = (row['pid'], row['tcov'], row['qcov'])
    metrics = {}
    for contig, representative in membership.items():
        if contig == representative:
            metrics[contig] = ('100', '100', '100')
        elif (representative, contig) in alignments:
            metrics[contig] = alignments[representative, contig]
        else:
            raise ValueError(f'Missing representative-to-contig ANI for {contig}')
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(FIELDS)
        for contig, length in lengths.items():
            representative = membership[contig]
            writer.writerow([contig, length, 'TRUE' if contig in selected else 'FALSE',
                             representative, lengths[representative], cluster_sizes[representative], *metrics[contig]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['fasta', 'clusters', 'representatives', 'ani', 'output']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    write_report(args.fasta, args.clusters, args.representatives, args.ani, args.output)


if __name__ == '__main__':
    main()
