#!/usr/bin/env python3
"""Write contig and summary reports for one sample."""
import argparse
import csv
from pathlib import Path

CONTIG_FIELDS = ['Viral contig', 'Sample', 'Assembler', 'Assembly reads', 'Percentage bacterial reads', 'All reads number', 'Reads used for assembly number', 'Total bp', 'Total corrected bp', 'Coverage', 'Coverage without bacteria', 'Flye circularity',
                 'Number of concatemers broken', 'DTR length removed']
ANNOTATION_FIELDS = ['CheckV gene count', 'CheckV viral genes', 'CheckV host genes',
                     'CheckV quality', 'MIUVIG quality', 'CheckV completeness', 'CheckV contamination',
                     'geNomad provirus', 'geNomad taxonomy']
CONTIG_FIELDS += ANNOTATION_FIELDS
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

def fasta_lengths(path):
    lengths = {}
    contig = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                contig = line[1:].split()[0]
                if contig in lengths:
                    raise ValueError(f'Duplicate corrected contig ID: {contig}')
                lengths[contig] = 0
            elif contig is None:
                raise ValueError(f'Sequence before FASTA header in {path}')
            else:
                lengths[contig] += len(''.join(line.split()))
    return lengths


def sample_report(sample, fasta, flye_info, concatemer_report, dtr_report, output, checkv=None, genomad=None, summary_output=None, assembler="flye", read_stats=None):
    read_fields = ['Assembly reads', 'Percentage bacterial reads', 'All reads number', 'Reads used for assembly number']
    read_values = dict.fromkeys(read_fields, 'NA')
    if read_stats:
        with open(read_stats) as handle:
            read_values = next(csv.DictReader(handle, delimiter='\t'))
    lengths = fasta_lengths(fasta)
    ids = set(lengths)
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
            if len(fields) < 5 or fields[0] in flye:
                raise ValueError(f'Invalid or duplicate assembly metadata in {flye_info}')
            flye[fields[0]] = fields
    rows = []
    for contig, corrected_length in lengths.items():
        if contig not in flye:
            raise ValueError(f'No assembly metadata for {contig}')
        metadata = flye[contig]
        original_length = int(metadata[1])
        circular = metadata[4]
        concatemer, dtr = concatemers[contig], dtrs[contig]
        if int(dtr['original_length']) != int(concatemer['corrected_length']) or int(dtr['corrected_length']) != corrected_length:
            raise ValueError(f'Correction lengths do not match for {contig}')
        if int(concatemer['original_length']) != original_length:
            raise ValueError(f'Assembly and correction lengths do not match for {contig}')
        if metadata[4] not in {'Y', 'N', 'NA'}:
            raise ValueError(f'Invalid assembly circularity for {contig}')
        float(metadata[2])  # Reject missing/malformed coverage rather than inventing a value.
        float(metadata[3])
        rows.append(dict(zip(CONTIG_FIELDS, [
            contig, sample, assembler, *[read_values[field] for field in read_fields], original_length, corrected_length, metadata[2], metadata[3], circular,
            int(concatemer['num_copies']) - 1 if concatemer['status'] == 'corrected' else 0,
            int(dtr['original_length']) - int(dtr['corrected_length']),
        ])))
    for row in rows:
        row.update(dict.fromkeys(ANNOTATION_FIELDS, 'NA'))
    if rows and (checkv or genomad):
        add_quality_metadata(rows, sample, checkv, genomad)
    write_tsv(output, CONTIG_FIELDS, rows)
    if summary_output is not None:
        write_tsv(summary_output, SAMPLE_FIELDS, [sample_summary(sample, rows)])


def sample_summary(sample, rows):
    return dict(zip(SAMPLE_FIELDS, [
        sample, len(rows), sum(int(row['Total bp']) for row in rows),
        sum(int(row['Total corrected bp']) for row in rows),
        'ready' if rows else 'no_viral_contigs',
    ]))


def indexed_tsv(path, key):
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    result = {row[key]: row for row in rows}
    if len(rows) != len(result):
        raise ValueError(f'Duplicate {key} in {path}')
    return result

def add_quality_metadata(rows, sample, checkv, genomad):
    if not all([checkv, genomad]):
        raise ValueError('Nonempty contig reports require CheckV and geNomad results')
    quality = indexed_tsv(checkv, 'contig_id')
    taxonomy = indexed_tsv(genomad, 'seq_name')
    columns = {'CheckV gene count': 'gene_count', 'CheckV viral genes': 'viral_genes',
               'CheckV host genes': 'host_genes', 'CheckV quality': 'checkv_quality',
               'MIUVIG quality': 'miuvig_quality', 'CheckV completeness': 'completeness',
               'CheckV contamination': 'contamination'}
    for row in rows:
        contig = row['Viral contig']
        original = contig.removeprefix(sample + '_')
        if contig not in quality or original not in taxonomy:
            raise ValueError(f'Missing CheckV or geNomad results for {contig}')
        row.update({name: quality[contig][source] for name, source in columns.items()})
        row['geNomad provirus'] = 'Yes' if taxonomy[original]['topology'] == 'Provirus' else 'No'
        row['geNomad taxonomy'] = taxonomy[original]['taxonomy']

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['sample', 'fasta', 'concatemer-report', 'dtr-report',
                 'checkv', 'genomad', 'output', 'summary-output']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument("--assembly-info", "--flye-info", dest="flye_info", required=True)
    parser.add_argument("--assembler", required=True)
    parser.add_argument("--read-stats", required=True)
    args = parser.parse_args()
    sample_report(args.sample, args.fasta, args.flye_info, args.concatemer_report,
                  args.dtr_report, args.output, args.checkv, args.genomad, args.summary_output, args.assembler, args.read_stats)


if __name__ == '__main__':
    main()
