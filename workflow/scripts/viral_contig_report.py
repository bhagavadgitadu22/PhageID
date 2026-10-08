#!/usr/bin/env python3
"""Report corrected viral contigs and select samples for downstream analysis."""
import argparse
import csv
from pathlib import Path
from Bio import SeqIO

CONTIG_FIELDS = ['Viral contig', 'Sample', 'Total bp', 'Total corrected bp', 'Coverage', 'Circular',
                 'Number of concatemers broken', 'DTR length removed']
ANNOTATION_FIELDS = ['CheckV gene count', 'CheckV viral genes', 'CheckV host genes',
                     'CheckV quality', 'MIUVIG quality', 'CheckV completeness', 'CheckV completeness_method',
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

def sample_report(sample, fasta, flye_info, concatemer_report, dtr_report, output, checkv=None, genomad=None):
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
        original_id = record.id[len(prefix):] if record.id.startswith(prefix) else ''
        flye_id, separator, region = original_id.partition('|provirus_')
        if flye_id not in flye:
            raise ValueError(f'No Flye metadata for {record.id}')
        metadata = flye[flye_id]
        original_length = int(metadata[1])
        circular = metadata[3]
        if separator:
            start, end = map(int, region.split('_'))
            if not 1 <= start <= end <= original_length:
                raise ValueError(f'Invalid provirus coordinates for {record.id}')
            original_length = end - start + 1
            circular = 'N'
        concatemer, dtr = concatemers[record.id], dtrs[record.id]
        if int(dtr['original_length']) != int(concatemer['corrected_length']) or int(dtr['corrected_length']) != len(record.seq):
            raise ValueError(f'Correction lengths do not match for {record.id}')
        if int(concatemer['original_length']) != original_length:
            raise ValueError(f'Flye and correction lengths do not match for {record.id}')
        if metadata[3] not in {'Y', 'N'}:
            raise ValueError(f'Invalid Flye circularity for {record.id}')
        float(metadata[2])  # Reject missing/malformed coverage rather than inventing a value.
        rows.append(dict(zip(CONTIG_FIELDS, [
            record.id, sample, original_length, len(record.seq), metadata[2], circular,
            int(concatemer['num_copies']) - 1 if concatemer['status'] == 'corrected' else 0,
            int(dtr['original_length']) - int(dtr['corrected_length']),
        ])))
    for row in rows:
        row.update(dict.fromkeys(ANNOTATION_FIELDS, 'NA'))
    if rows and (checkv or genomad):
        add_quality_metadata(rows, sample, checkv, genomad)
    write_tsv(output, CONTIG_FIELDS, rows)

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
               'CheckV completeness_method': 'completeness_method'}
    for row in rows:
        contig = row['Viral contig']
        original = contig.removeprefix(sample + '_')
        if contig not in quality or original not in taxonomy:
            raise ValueError(f'Missing CheckV or geNomad results for {contig}')
        row.update({name: quality[contig][source] for name, source in columns.items()})
        row['geNomad provirus'] = 'Yes' if taxonomy[original]['topology'] == 'Provirus' else 'No'
        row['geNomad taxonomy'] = taxonomy[original]['taxonomy']

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

def build_reports(samples, results_dir, output_dir):
    reports = []
    for sample, fasta, flye, concatemer, dtr, checkv, genomad in samples:
        output = Path(results_dir) / sample / 'assembly_stats.tsv'
        sample_report(sample, fasta, flye, concatemer, dtr, output, checkv, genomad)
        reports.append((sample, output))
    combine_reports(reports, output_dir)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-input', action='append', nargs=7, default=[],
                        metavar=('SAMPLE', 'FASTA', 'FLYE', 'CONCATEMER', 'DTR', 'CHECKV', 'GENOMAD'))
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    build_reports(args.sample_input, args.results_dir, args.output_dir)

if __name__ == '__main__':
    main()
