#!/usr/bin/env python3
"""Report corrected viral contigs and select samples for downstream analysis."""
import argparse
import csv
from pathlib import Path
from Bio import SeqIO

CONTIG_FIELDS = ['Viral contig', 'Sample', 'Total bp', 'Total corrected bp', 'Coverage', 'Circular',
                 'Number of concatemers broken', 'DTR length removed']
SAMPLE_FIELDS = ['Sample', 'Viral contigs', 'Total bp', 'Total corrected bp', 'Status']


def read_corrections(path):
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    result = {row['contig_id']: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f'Duplicate correction IDs in {path}')
    return result


def write_tsv(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sample_report(sample, fasta, flye_info, concatemer_report, dtr_report, output):
    records = list(SeqIO.parse(fasta, 'fasta'))
    ids = {record.id for record in records}
    if len(ids) != len(records):
        raise ValueError(f'Duplicate corrected contig IDs for {sample}')
    concatemers = read_corrections(concatemer_report)
    dtrs = read_corrections(dtr_report)
    if ids != set(concatemers) or ids != set(dtrs):
        raise ValueError(f'Correction report IDs do not match corrected FASTA for {sample}')
    flye = {}
    with open(flye_info) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4 or fields[0] in flye:
                raise ValueError(f'Invalid or duplicate Flye metadata in {flye_info}')
            flye[fields[0]] = fields
    rows = []
    prefix = sample + '_'
    for record in records:
        # rename_contigs prefixes the original Flye ID with the full sample ID.
        if not record.id.startswith(prefix) or record.id[len(prefix):] not in flye:
            raise ValueError(f'No Flye metadata for {record.id}')
        metadata = flye[record.id[len(prefix):]]
        concatemer, dtr = concatemers[record.id], dtrs[record.id]
        if int(dtr['original_length']) != int(concatemer['corrected_length']) or int(dtr['corrected_length']) != len(record.seq):
            raise ValueError(f'Correction lengths do not match for {record.id}')
        if int(concatemer['original_length']) != int(metadata[1]):
            raise ValueError(f'Flye and correction lengths do not match for {record.id}')
        if metadata[3] not in {'Y', 'N'}:
            raise ValueError(f'Invalid Flye circularity for {record.id}')
        float(metadata[2])  # Reject missing/malformed coverage rather than inventing a value.
        rows.append(dict(zip(CONTIG_FIELDS, [
            record.id, sample, int(metadata[1]), len(record.seq), metadata[2], metadata[3],
            int(concatemer['num_copies']) - 1 if concatemer['status'] == 'corrected' else 0,
            int(dtr['original_length']) - int(dtr['corrected_length']),
        ])))
    write_tsv(output, CONTIG_FIELDS, rows)


def combine_reports(sample_reports, output_dir):
    all_rows, summaries = [], []
    seen = set()
    for sample, report in sample_reports:
        with open(report, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if reader.fieldnames != CONTIG_FIELDS:
                raise ValueError(f'Unexpected contig report columns: {report}')
            rows = list(reader)
        for row in rows:
            if row['Sample'] != sample or row['Viral contig'] in seen:
                raise ValueError(f'Duplicate contig or incorrect sample in {report}')
            seen.add(row['Viral contig'])
        all_rows.extend(rows)
        summaries.append(dict(zip(SAMPLE_FIELDS, [
            sample, len(rows), sum(int(row['Total bp']) for row in rows),
            sum(int(row['Total corrected bp']) for row in rows),
            'ready' if rows else 'no_viral_contigs',
        ])))
    output = Path(output_dir)
    write_tsv(output/'viral_contigs.tsv', CONTIG_FIELDS, all_rows)
    write_tsv(output/'samples.tsv', SAMPLE_FIELDS, summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    sample = commands.add_parser('sample')
    for name in ['sample', 'fasta', 'flye-info', 'concatemer-report', 'dtr-report', 'output']:
        sample.add_argument('--' + name, required=True)
    combine = commands.add_parser('combine')
    combine.add_argument('--sample-report', action='append', nargs=2, default=[])
    combine.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    if args.command == 'sample':
        sample_report(args.sample, args.fasta, args.flye_info, args.concatemer_report, args.dtr_report, args.output)
    else:
        combine_reports(args.sample_report, args.output_dir)


if __name__ == '__main__':
    main()
