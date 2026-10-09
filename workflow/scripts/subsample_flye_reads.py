#!/usr/bin/env python3
"""Reproducibly subsample whole FASTQ records to a target number of bases."""
import argparse
import csv
import random


def fastq_records(path):
    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence, separator, quality = [handle.readline() for _ in range(3)]
            length = len(sequence.rstrip('\r\n'))
            if (not header.startswith('@') or not separator.startswith('+') or
                    not length or len(quality.rstrip('\r\n')) != length):
                raise ValueError(f'Invalid or truncated FASTQ record in {path}')
            yield (header, sequence, separator, quality), length


def subsample(reads, output, report, genome_size=100000, coverage=1000, seed=42):
    if genome_size <= 0 or coverage <= 0:
        raise ValueError('Genome size and target coverage must be positive')
    lengths = [length for _, length in fastq_records(reads)]
    if not lengths:
        raise ValueError('No reads available for Flye')
    target = genome_size * coverage
    input_bases = sum(lengths)
    selected = None
    if input_bases > target:
        indices = list(range(len(lengths)))
        random.Random(seed).shuffle(indices)
        selected = set()
        bases = 0
        for index in indices:
            selected.add(index)
            bases += lengths[index]
            if bases >= target:
                break
    output_bases = output_reads = 0
    with open(output, 'w') as handle:
        for index, (record, length) in enumerate(fastq_records(reads)):
            if selected is None or index in selected:
                handle.writelines(record)
                output_bases += length
                output_reads += 1
    with open(report, 'w', newline='') as handle:
        fields = ['Genome size estimate bp', 'Target coverage', 'Target bp', 'Input reads',
                  'Input bp', 'Selected reads', 'Selected bp', 'Estimated coverage', 'Seed']
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerow([genome_size, coverage, target, len(lengths), input_bases,
                         output_reads, output_bases, output_bases / genome_size, seed])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['reads', 'output', 'report']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--genome-size', type=int, default=100000)
    parser.add_argument('--coverage', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    subsample(args.reads, args.output, args.report, args.genome_size, args.coverage, args.seed)


if __name__ == '__main__':
    main()
