#!/usr/bin/env python3
"""Prepare matching proteins and gene-to-genome mappings for dereplicated viruses."""
import argparse
import csv
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO

def prepare(representatives, genbanks, proteins, gene_map, genome_lengths):
    records = list(SeqIO.parse(representatives, 'fasta'))
    lengths = {record.id: len(record.seq) for record in records}
    if len(records) != len(lengths):
        raise ValueError('Duplicate representative IDs')
    sequences, mappings, seen = [], [], set()
    for path in genbanks:
        for record in SeqIO.parse(path, 'genbank'):
            if record.id not in lengths:
                continue
            if record.id in seen:
                raise ValueError(f'Duplicate representative annotation: {record.id}')
            seen.add(record.id)
            if len(record.seq) != lengths[record.id]:
                raise ValueError(f'Annotation length differs for {record.id}')
            count = 0
            for feature in record.features:
                if feature.type != 'CDS':
                    continue
                translation = feature.qualifiers.get('translation', [])
                if not translation or not translation[0].strip():
                    raise ValueError(f'Missing CDS translation for {record.id}')
                count += 1
                protein = quote(record.id, safe='_-.') + '__CDS_' + str(count)
                sequences.append((protein, ''.join(translation[0].split())))
                mappings.append((protein, record.id))
            if count == 0:
                raise ValueError(f'No annotated proteins for {record.id}')
    if seen != set(lengths):
        raise ValueError(f'Missing representative annotations: {sorted(set(lengths) - seen)}')
    for path in [proteins, gene_map, genome_lengths]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(proteins, 'w') as handle:
        for protein, sequence in sequences:
            handle.write(f'>{protein}\n{sequence}\n')
    with open(gene_map, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['protein_id', 'genome_id'])
        writer.writerows(mappings)
    with open(genome_lengths, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['genome_id', 'length'])
        writer.writerows(lengths.items())

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genbanks', nargs='*', default=[])
    for name in ['representatives', 'proteins', 'gene-map', 'genome-lengths']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    prepare(args.representatives, args.genbanks, args.proteins, args.gene_map, args.genome_lengths)

if __name__ == '__main__':
    main()
