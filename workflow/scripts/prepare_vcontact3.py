#!/usr/bin/env python3
"""Prepare matching proteins and gene-to-genome mappings for dereplicated viruses."""
import argparse
from glob import glob
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO
import pandas as pd


def resolve_genbanks(inputs):
    """Expand explicit files, directories, and quoted shell globs."""
    paths = []
    seen = set()
    for value in inputs:
        pattern = str(Path(value).expanduser())
        matches = sorted(glob(pattern, recursive=True))
        if not matches:
            raise ValueError(f'No GenBank files match: {value}')
        for match in matches:
            path = Path(match)
            candidates = (sorted(p for p in path.rglob('*')
                                 if p.is_file() and p.suffix.lower() in {'.gbk', '.gbff', '.gb'})
                          if path.is_dir() else [path])
            if not candidates:
                raise ValueError(f'No GenBank files found in: {path}')
            for candidate in candidates:
                if not candidate.is_file():
                    raise ValueError(f'Not a GenBank file: {candidate}')
                canonical = candidate.resolve()
                if canonical not in seen:
                    paths.append(candidate)
                    seen.add(canonical)
    return paths


def prepare(representatives, genbanks, proteins, gene_map, genome_lengths):
    records = list(SeqIO.parse(representatives, 'fasta'))
    representatives_by_id = {record.id: record for record in records}
    if len(records) != len(representatives_by_id):
        raise ValueError('Duplicate representative IDs')
    lengths = {record.id: len(record.seq) for record in records}
    genbank_paths = resolve_genbanks(genbanks)
    if records and not genbank_paths:
        raise ValueError('No GenBank files supplied for nonempty representatives')
    sequences, mappings, seen = [], [], set()
    for path in genbank_paths:
        with open(path) as handle:
            for record in SeqIO.parse(handle, 'genbank'):
                if record.id not in lengths:
                    continue
                if record.id in seen:
                    raise ValueError(f'Duplicate representative annotation: {record.id} (in {path})')
                seen.add(record.id)
                if len(record.seq) != lengths[record.id]:
                    raise ValueError(f'Annotation length differs for {record.id} in {path}')
                if record.seq.upper() != representatives_by_id[record.id].seq.upper():
                    raise ValueError(f'Annotation sequence differs for {record.id} in {path}')
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
                    product = feature.qualifiers.get('product', ['None'])
                    keyword = product[0].strip() if product and product[0].strip() else 'None'
                    mappings.append((protein, keyword, record.id))
                if count == 0:
                    raise ValueError(f'No annotated proteins for {record.id}')
    if seen != set(lengths):
        raise ValueError(f'Missing representative annotations: {sorted(set(lengths) - seen)} '
                         f'(searched {len(genbank_paths)} GenBank files: '
                         f'{[str(path) for path in genbank_paths]})')
    for path in [proteins, gene_map, genome_lengths]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(proteins, 'w') as handle:
        for protein, sequence in sequences:
            handle.write(f'>{protein}\n{sequence}\n')
    pd.DataFrame(mappings, columns=['protein_id', 'keywords', 'genome_id']).to_parquet(gene_map, index=False)
    pd.DataFrame(lengths.items(), columns=['genome_id', 'length']).to_parquet(genome_lengths, index=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genbanks', nargs='*', default=[],
                        help='GenBank files, directories, or glob patterns (quoted patterns are accepted)')
    for name in ['representatives', 'proteins', 'gene-map', 'genome-lengths']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    prepare(args.representatives, args.genbanks, args.proteins, args.gene_map, args.genome_lengths)

if __name__ == '__main__':
    main()
