#!/usr/bin/env python3
"""Persist a resolved assembly choice and validate sources before publication."""
import argparse
import csv
import json
import shutil
import tempfile
from pathlib import Path

from assembly_read_stats import selected_stats
from subsample_flye_reads import fastq_records
from viral_contig_report import fasta_lengths

SOURCE_FILES = {
    'corrected': 'circularisation/circular_viruses.fasta',
    'concatemer_report': 'circularisation/breaking_concatemers_report.csv',
    'dtr_report': 'circularisation/breaking_dtr_report.csv',
    'checkv_quality': 'checkv/quality_summary.tsv',
    'fasta': 'genomad/geNomad_assembly/assembly_summary/assembly_virus.fna',
    'summary': 'genomad/geNomad_assembly/assembly_summary/assembly_virus_summary.tsv',
    'assembly': 'assembly.fasta',
    'status': 'assembly_status.txt',
}
DESTINATION_KEYS = set(SOURCE_FILES) | {'flye_info', 'assembler', 'read_stats'}


def table(path, columns, delimiter='\t'):
    with open(path, newline='') as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not set(columns) <= set(reader.fieldnames or []):
            raise ValueError(f'Invalid table header: {path}; expected {columns}')
        return list(reader)


def validate(manifest):
    sample, assembler, readset = (manifest[key] for key in ['sample', 'assembler', 'readset'])
    if assembler not in {'flye', 'spades', 'autocycler'} or readset not in {'filtered', 'all'}:
        raise ValueError('Invalid assembler or read set in selection manifest')
    sources = {key: Path(value) for key, value in manifest['sources'].items()}
    expected_keys = set(SOURCE_FILES) | {'used_reads', 'filter_report'}
    if assembler == 'flye':
        expected_keys.add('flye_info')
    if set(sources) != expected_keys:
        raise ValueError('Missing or unexpected source keys in selection manifest')
    candidate = sources['assembly'].parent
    if (candidate.name, candidate.parent.name, candidate.parent.parent.name,
            candidate.parent.parent.parent.name) != (assembler, readset, 'assembly_attempts', sample):
        raise ValueError('Assembly path does not match the selected sample/assembler/read set')
    for key, relative in SOURCE_FILES.items():
        if sources[key] != candidate / relative:
            raise ValueError(f'Unexpected {key} source path: {sources[key]}')
    sample_dir = candidate.parent.parent.parent
    if sources['filter_report'] != sample_dir / 'reads/host_filtering_stats.tsv':
        raise ValueError('Unexpected filtering report path')
    if assembler == 'flye' and sources['flye_info'] != candidate / 'assembly_info.txt':
        raise ValueError('Unexpected Flye metadata path')
    expected_reads = (sample_dir / 'reads' / f'flye.{sample}.{readset}.fastq' if assembler == 'flye'
                      else sample_dir / 'reads' / f'cleaned.{sample}.fastq' if readset == 'filtered'
                      else None)
    if expected_reads is not None and sources['used_reads'] != expected_reads:
        raise ValueError('Unexpected assembly read path')
    if expected_reads is None and sources['used_reads'] not in {
            sample_dir / 'reads' / f'porechop.{sample}.fastq',
            sample_dir / 'reads' / f'fastp.{sample}.fastq'}:
        raise ValueError('Unexpected all-read input path')
    for path in sources.values():
        if not path.is_file():
            raise ValueError(f'Missing selection source: {path}')
    status = sources['status'].read_text().strip()
    if status not in {'success', 'failed'}:
        raise ValueError(f'Invalid assembly status: {sources["status"]}')
    lengths = {key: fasta_lengths(sources[key]) for key in ['assembly', 'corrected', 'fasta']}
    if any(length <= 0 for records in lengths.values() for length in records.values()):
        raise ValueError('Empty sequence record in selected FASTA')
    if status == 'failed' and any(lengths.values()):
        raise ValueError('Failed assembly cannot contain sequences')
    if status == 'success' and not lengths['assembly']:
        raise ValueError('Successful assembly has no sequences')
    ids = set(lengths['corrected'])
    raw_ids = set(lengths['fasta'])
    original_ids = {contig.removeprefix(sample + '_') for contig in ids}
    if not (ids <= raw_ids or original_ids <= raw_ids):
        raise ValueError('Corrected contigs are missing from the geNomad FASTA')
    for key in ['concatemer_report', 'dtr_report']:
        rows = table(sources[key], ['contig_id', 'original_length', 'corrected_length'], ',')
        if {row['contig_id'] for row in rows} != ids or len(rows) != len(ids):
            raise ValueError(f'Correction IDs do not match corrected FASTA: {sources[key]}')
    rows = table(sources['checkv_quality'], ['contig_id', 'checkv_quality'])
    if {row['contig_id'] for row in rows} != ids or len(rows) != len(ids):
        raise ValueError('CheckV IDs do not match corrected FASTA')
    table(sources['summary'], ['seq_name', 'topology', 'taxonomy'])
    rows = table(sources['filter_report'], ['All reads number', 'Bacterial free reads', 'Percentage bacterial reads'])
    if len(rows) != 1:
        raise ValueError('Expected one filtering report row')
    if assembler == 'flye' and status == 'success':
        table(sources['flye_info'], ['#seq_name', 'length', 'cov.', 'circ.'])
    if next(fastq_records(sources['used_reads']), None) is None:
        raise ValueError('Selected assembly reads are empty')
    return sources


def create(sample, inputs, output):
    paths = [Path(path).resolve() for path in inputs]
    assemblies = [path for path in paths if path.name == 'assembly.fasta']
    reads = [path for path in paths if path.suffix == '.fastq']
    if len(assemblies) != 1 or len(reads) != 1:
        raise ValueError('Selection inputs are unresolved or ambiguous; expected one assembly FASTA and one read file')
    candidate = assemblies[0].parent
    sources = {key: candidate / relative for key, relative in SOURCE_FILES.items()}
    sources['used_reads'] = reads[0]
    sources['filter_report'] = candidate.parent.parent.parent / 'reads/host_filtering_stats.tsv'
    if candidate.name == 'flye':
        sources['flye_info'] = candidate / 'assembly_info.txt'
    if not set(sources.values()) <= set(paths):
        raise ValueError('Selection sources were not all declared as completed inputs')
    manifest = {'sample': sample, 'assembler': candidate.name, 'readset': candidate.parent.name,
                'sources': {key: str(path) for key, path in sources.items()}}
    validate(manifest)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + '\n')


def publish(sample, manifest_path, destinations):
    manifest = json.loads(Path(manifest_path).read_text())
    if manifest['sample'] != sample:
        raise ValueError('Selection manifest belongs to another sample')
    if set(destinations) != DESTINATION_KEYS:
        raise ValueError('Missing or unexpected publication destinations')
    sources = validate(manifest)
    # Count and validate all read records before touching any published outputs.
    with tempfile.TemporaryDirectory(prefix='phageid_selection_') as folder:
        read_stats = Path(folder) / 'read_stats.tsv'
        selected_stats(sources['filter_report'], sources['used_reads'], manifest['readset'], read_stats)
        for path in destinations.values():
            Path(path).parent.mkdir(parents=True, exist_ok=True)
        for key in SOURCE_FILES:
            shutil.copyfile(sources[key], destinations[key])
        if manifest['assembler'] == 'flye':
            shutil.copyfile(sources['flye_info'], destinations['flye_info'])
        else:
            Path(destinations['flye_info']).write_text('#seq_name\tlength\tcov.\tcirc.\n')
        Path(destinations['assembler']).write_text(manifest['assembler'] + '\n')
        shutil.copyfile(read_stats, destinations['read_stats'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode', required=True)
    create_parser = modes.add_parser('create')
    create_parser.add_argument('--sample', required=True)
    create_parser.add_argument('--inputs', nargs='+', required=True)
    create_parser.add_argument('--output', required=True)
    publish_parser = modes.add_parser('publish')
    publish_parser.add_argument('--sample', required=True)
    publish_parser.add_argument('--manifest', required=True)
    publish_parser.add_argument('--destination', nargs=2, action='append', required=True)
    args = parser.parse_args()
    if args.mode == 'create':
        create(args.sample, args.inputs, args.output)
    else:
        destinations = dict(args.destination)
        if len(destinations) != len(args.destination):
            parser.error('Duplicate destination keys')
        publish(args.sample, args.manifest, destinations)


if __name__ == '__main__':
    main()
