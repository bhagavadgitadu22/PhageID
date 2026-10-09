#!/usr/bin/env python3
"""Count preprocessed, host-filtered and actual assembly input reads."""
import argparse
import csv
from subsample_flye_reads import fastq_records

FIELDS = ['Assembly reads', 'Percentage bacterial reads', 'All reads number', 'Reads used for assembly number']

def count_reads(path):
    return sum(1 for _ in fastq_records(path))

def write_row(output, row):
    with open(output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)

def filtering_stats(all_reads, filtered_reads, has_host, output):
    total, filtered = count_reads(all_reads), count_reads(filtered_reads)
    if total == 0 or filtered > total:
        raise ValueError('Invalid read counts after host filtering')
    percentage = round(100 * (total - filtered) / total, 2) if has_host else 'NA'
    write_row(output, {'All reads number': total, 'Bacterial free reads': filtered,
                      'Percentage bacterial reads': percentage})

def selected_stats(filter_report, used_reads, readset, output):
    with open(filter_report) as handle:
        stats = next(csv.DictReader(handle, delimiter='\t'))
    used = count_reads(used_reads)
    available = int(stats['Bacterial free reads'] if readset == 'filtered' else stats['All reads number'])
    if not 0 < used <= available:
        raise ValueError('Invalid number of assembly input reads')
    write_row(output, dict(zip(FIELDS, ['bacterial free' if readset == 'filtered' else 'all',
              stats['Percentage bacterial reads'], stats['All reads number'], used])))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode', required=True)
    filtering = modes.add_parser('filter')
    for name in ['all-reads', 'filtered-reads', 'output']:
        filtering.add_argument('--' + name, required=True)
    filtering.add_argument('--has-host', type=int, choices=[0, 1], required=True)
    selected = modes.add_parser('selected')
    for name in ['filter-report', 'used-reads', 'output']:
        selected.add_argument('--' + name, required=True)
    selected.add_argument('--readset', choices=['filtered', 'all'], required=True)
    args = parser.parse_args()
    if args.mode == 'filter':
        filtering_stats(args.all_reads, args.filtered_reads, args.has_host, args.output)
    else:
        selected_stats(args.filter_report, args.used_reads, args.readset, args.output)

if __name__ == '__main__':
    main()
