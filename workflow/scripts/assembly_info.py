#!/usr/bin/env python3
"""Report circular-BAM coverage on corrected viruses and Flye assembly metadata."""
import argparse
import csv
from viral_contig_report import fasta_lengths, read_corrections


def bam_depths(path, lengths):
    import pysam
    totals = dict.fromkeys(lengths, 0)
    with pysam.AlignmentFile(path, 'rb', check_sq=False) as bam:
        if dict(zip(bam.references, bam.lengths)) != lengths:
            raise ValueError(f'BAM references do not match corrected contigs: {path}')
        circular = 'theBIGbam:circular=true' in bam.header.to_dict().get('CO', [])
        for read in bam.fetch(until_eof=True):
            # Same basic exclusion flags as samtools coverage; include supplementary.
            if read.flag & (4 | 256 | 512 | 1024):
                continue
            name = read.reference_name
            length = lengths[name]
            position = read.reference_start
            for operation, count in read.cigartuples or []:
                if operation in (0, 7, 8):  # M, =, X contribute aligned bases.
                    totals[name] += count if circular else max(0, min(position + count, length) - max(position, 0))
                    position += count
                elif operation in (2, 3):  # D and N consume reference but add no depth.
                    position += count
    return {name: totals[name] / length for name, length in lengths.items()}


def write_info(sample, fasta, corrected, concatemer_report, bam, filtered_bam, flye_info, output):
    originals = fasta_lengths(fasta)
    lengths = fasta_lengths(corrected)
    corrections = read_corrections(concatemer_report)
    if set(corrections) != set(lengths):
        raise ValueError('Correction records do not match corrected contigs')
    depths = bam_depths(bam, lengths)
    filtered_depths = bam_depths(filtered_bam, lengths)
    with open(flye_info, newline='') as handle:
        metadata = {row['#seq_name']: row for row in csv.DictReader(handle, delimiter='\t')}
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['#seq_name', 'length', 'cov.', 'cov_without_bacteria', 'flye_circularity'])
        for name in lengths:
            original = name.removeprefix(sample + '_')
            parent, separator, region = original.partition('|provirus_')
            original_length = originals[parent]
            if separator:
                start, end = map(int, region.split('_'))
                if not 1 <= start <= end <= original_length:
                    raise ValueError(f'Invalid provirus coordinates: {name}')
                original_length = end - start + 1
            if int(corrections[name]['original_length']) != original_length:
                raise ValueError(f'Correction length does not match original assembly: {name}')
            circular = metadata[parent]['circ.'] if parent in metadata and not separator else 'NA'
            if circular not in {'Y', 'N', 'NA'}:
                raise ValueError(f'Invalid Flye circularity: {parent}')
            writer.writerow([name, original_length, depths[name], filtered_depths[name], circular])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['sample', 'fasta', 'corrected', 'concatemer-report', 'bam', 'filtered-bam', 'flye-info', 'output']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    write_info(args.sample, args.fasta, args.corrected, args.concatemer_report,
               args.bam, args.filtered_bam, args.flye_info, args.output)


if __name__ == '__main__':
    main()
