#!/usr/bin/env python3
"""Write contig and summary reports for one sample."""
import argparse
import csv
from pathlib import Path

CONTIG_FIELDS = ['Viral_contig', 'Sample', 'Number_reads', 'Number_reads_without_host', 'Assembler', 'Assembly_method', 'Assembly_used_reads', 'Total_bp', 'Total_corrected_bp', 'Coverage', 'Coverage_without_bacteria', 'Flye_circularity',
                 'Number_of_concatemers_broken', 'DTR_length_removed', 'Manual_correction', 'Manual_length_change_bp']
ANNOTATION_FIELDS = ['CheckV_gene_count', 'CheckV_viral_genes', 'CheckV_host_genes',
                     'CheckV_quality', 'MIUVIG_quality', 'CheckV_completeness', 'CheckV_contamination',
                     'geNomad_provirus', 'geNomad_taxonomy']
CONTIG_FIELDS += ANNOTATION_FIELDS
SAMPLE_FIELDS = ['Sample', 'Number_reads', 'Number_reads_without_host', 'Assembler', 'Assembly_method', 'Assembly_used_reads', 'Viral_contigs', 'Total_bp', 'Total_corrected_bp', 'Status']

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


def sample_report(sample, fasta, flye_info, concatemer_report, dtr_report, output, checkv=None, genomad=None, summary_output=None, assembler="flye", read_stats=None, assembly_status="success", filter_report=None, manual_report=None):
    read_fields = ['Assembly reads', 'Percentage bacterial reads', 'All reads number', 'Reads used for assembly number']
    read_values = dict.fromkeys(read_fields, 'NA')
    if read_stats:
        with open(read_stats) as handle:
            read_values = next(csv.DictReader(handle, delimiter='\t'))
    filtering = {}
    if filter_report:
        with open(filter_report) as handle:
            filtering = next(csv.DictReader(handle, delimiter='\t'))
    assembly_values = [read_values['All reads number'], filtering.get('Bacterial free reads', read_values.get('Bacterial free reads', 'NA')), assembler, read_values['Assembly reads'], read_values['Reads used for assembly number']]
    lengths = fasta_lengths(fasta)
    ids = set(lengths)
    concatemers = read_corrections(concatemer_report)
    dtrs = read_corrections(dtr_report)
    if ids != set(concatemers) or ids != set(dtrs):
        raise ValueError(f'Correction report IDs do not match corrected FASTA for {sample}')
    manual = indexed_tsv(manual_report, 'Viral_contig') if manual_report else {}
    if manual_report and set(manual) != ids:
        raise ValueError(f'Manual correction report IDs do not match corrected FASTA for {sample}')
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
        auto_length = int(dtr['corrected_length'])
        if int(dtr['original_length']) != int(concatemer['corrected_length']):
            raise ValueError(f'Correction lengths do not match for {contig}')
        manual_change = manual.get(contig)
        if manual_change:
            if (int(manual_change['Auto_corrected_bp']) != auto_length
                    or int(manual_change['Final_bp']) != corrected_length
                    or int(manual_change['Manual_length_change_bp']) != corrected_length - auto_length
                    or manual_change['Manual_correction'] not in {'Yes', 'No'}):
                raise ValueError(f'Manual correction lengths do not match for {contig}')
            if manual_change['Manual_correction'] == 'No' and corrected_length != auto_length:
                raise ValueError(f'Manual correction status conflicts with length for {contig}')
        elif auto_length != corrected_length:
            raise ValueError(f'Correction lengths do not match for {contig}')
        if int(concatemer['original_length']) != original_length:
            raise ValueError(f'Assembly and correction lengths do not match for {contig}')
        if metadata[4] not in {'Y', 'N', 'NA'}:
            raise ValueError(f'Invalid assembly circularity for {contig}')
        float(metadata[2])  # Reject missing/malformed coverage rather than inventing a value.
        float(metadata[3])
        rows.append(dict(zip(CONTIG_FIELDS, [
            contig, sample, *assembly_values, original_length, corrected_length, metadata[2], metadata[3], circular,
            int(concatemer['num_copies']) - 1 if concatemer['status'] == 'corrected' else 0,
            int(dtr['original_length']) - auto_length,
            manual_change['Manual_correction'] if manual_change else 'No',
            corrected_length - auto_length,
        ])))
    for row in rows:
        row.update(dict.fromkeys(ANNOTATION_FIELDS, 'NA'))
    if rows and (checkv or genomad):
        add_quality_metadata(rows, sample, checkv, genomad)
    write_tsv(output, CONTIG_FIELDS, rows)
    if summary_output is not None:
        write_tsv(summary_output, SAMPLE_FIELDS, [sample_summary(sample, rows, assembly_status, assembly_values)])


def sample_summary(sample, rows, assembly_status="success", assembly_values=None):
    if assembly_status not in {"success", "failed"}:
        raise ValueError(f'Invalid assembly status: {assembly_status}')
    if assembly_status == "failed" and rows:
        raise ValueError('Failed assembly cannot have viral contigs')
    if assembly_values is None:
        assembly_values = [rows[0].get(field, 'NA') if rows else 'NA' for field in CONTIG_FIELDS[2:7]]
    return dict(zip(SAMPLE_FIELDS, [
        sample, *assembly_values, len(rows), sum(int(row['Total_bp']) for row in rows),
        sum(int(row['Total_corrected_bp']) for row in rows),
        'ready' if rows else 'assembly_failed' if assembly_status == 'failed' else 'no_viral_contigs',
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
    columns = {'CheckV_gene_count': 'gene_count', 'CheckV_viral_genes': 'viral_genes',
               'CheckV_host_genes': 'host_genes', 'CheckV_quality': 'checkv_quality',
               'MIUVIG_quality': 'miuvig_quality', 'CheckV_completeness': 'completeness',
               'CheckV_contamination': 'contamination'}
    for row in rows:
        contig = row['Viral_contig']
        original = contig.removeprefix(sample + '_')
        if contig not in quality or original not in taxonomy:
            raise ValueError(f'Missing CheckV or geNomad results for {contig}')
        row.update({name: quality[contig][source] for name, source in columns.items()})
        row['geNomad_provirus'] = 'Yes' if taxonomy[original]['topology'] == 'Provirus' else 'No'
        row['geNomad_taxonomy'] = taxonomy[original]['taxonomy']

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['sample', 'fasta', 'concatemer-report', 'dtr-report',
                 'checkv', 'genomad', 'output', 'summary-output']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument("--assembly-info", "--flye-info", dest="flye_info", required=True)
    parser.add_argument("--assembler", required=True)
    parser.add_argument("--read-stats", required=True)
    parser.add_argument("--filter-report", required=True)
    parser.add_argument("--assembly-status", required=True)
    parser.add_argument("--manual-report")
    args = parser.parse_args()
    assembly_status = Path(args.assembly_status).read_text().strip()
    sample_report(args.sample, args.fasta, args.flye_info, args.concatemer_report,
                  args.dtr_report, args.output, args.checkv, args.genomad, args.summary_output, args.assembler, args.read_stats, assembly_status, args.filter_report, args.manual_report)


if __name__ == '__main__':
    main()
