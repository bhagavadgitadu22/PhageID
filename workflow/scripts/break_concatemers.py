#!/usr/bin/env python
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description="Detect and correct concatemers from BLAST self-alignment")
    parser.add_argument('--blast', type=str, required=True, help="Path to input BLAST file")
    parser.add_argument('--fasta', type=str, required=True, help="Path to input FASTA file")
    parser.add_argument('--out_report', type=str, required=True, help="Path to output CSV report file")
    parser.add_argument('--out_fasta', type=str, required=True, help="Path to output corrected FASTA file")
    parser.add_argument('--min_identity', type=float, default=80, help="Minimum percent identity (default: 80)")
    parser.add_argument('--min_repeat', type=float, default=80, help="Repeats are considered when they represent at least x% of the biggest repeat (default: 80)")
    parser.add_argument('--min_coverage', type=float, default=80, help="Minimum coverage percent (default: 80)")
    return parser.parse_args()

args = parse_arguments()

# Parse BLAST results - only one sequence expected
# First pass: collect all hits above identity threshold
all_hits = []
with open(args.blast) as f:
    for line in f:
        fields = line.strip().split('\t')
        qseqid, sseqid, pident, length, qstart, qend, qlen, sstart, send, slen = fields
        
        # Skip self-hits at same position
        if qstart == sstart:
            continue
        
        pident = float(pident)
        length = int(length)
        if pident >= args.min_identity:
            all_hits.append({
                'qstart': int(qstart),
                'qend': int(qend),
                'sstart': int(sstart),
                'send': int(send),
                'length': length,
                'pident': pident,
                'qlen': int(qlen),
                'orientation': 'forward' if int(sstart) < int(send) else 'reverse'
            })

# Second pass: filter by min_repeat threshold relative to biggest repeat
blast_hits = []
if all_hits:
    max_length = max(hit['length'] for hit in all_hits)
    min_length_threshold = max_length * (args.min_repeat / 100.0)
    blast_hits = [hit for hit in all_hits if hit['length'] >= min_length_threshold]

# Process the single sequence
record = next(SeqIO.parse(args.fasta, "fasta"))
seq_len = len(record.seq)

# Result dictionary for CSV output
result = {
    'contig_id': record.id,
    'original_length': seq_len,
    'num_hits': 0,
    'total_repeat_bp': 0,
    'coverage_ratio': 0.0,
    'forward_hits': 0,
    'reverse_hits': 0,
    'concatemer_detected': False,
    'repeat_unit_size': 0,
    'num_copies': 0,
    'corrected_length': seq_len,
    'status': 'no_repeats'
}

if not blast_hits:
    result['status'] = 'no_repeats'
else:
    # Sum total bp in repeats
    total_repeat_bp = sum(hit['length'] for hit in blast_hits)
    coverage_ratio = total_repeat_bp / seq_len
    
    result['num_hits'] = len(blast_hits)
    result['total_repeat_bp'] = total_repeat_bp
    result['coverage_ratio'] = coverage_ratio
    
    # Check for inverted repeats (reverse complement)
    reverse_hits = [h for h in blast_hits if h['orientation'] == 'reverse']
    forward_hits = [h for h in blast_hits if h['orientation'] == 'forward']
    
    result['forward_hits'] = len(forward_hits)
    result['reverse_hits'] = len(reverse_hits)
    
    if coverage_ratio < args.min_coverage / 100.0:
        result['status'] = 'below_threshold'
    else:
        # Check for direct repeats (concatemers)
        if forward_hits:
            # After filtering, all repeats are similar size
            # Use the biggest one as the repeat unit size
            repeat_unit = max(hit['length'] for hit in forward_hits)
            num_copies = round(seq_len / repeat_unit)
            
            if num_copies >= 2:
                result['concatemer_detected'] = True
                result['repeat_unit_size'] = repeat_unit
                result['num_copies'] = num_copies
                result['corrected_length'] = repeat_unit
                result['status'] = 'corrected'
                
                # Extract single copy
                record.seq = record.seq[:repeat_unit]
                record.description = f"{record.description} | corrected from {num_copies}x concatemer"
            else:
                result['status'] = 'repeats_unclear'
        
        # Check for inverted repeats
        elif reverse_hits and not forward_hits:
            result['status'] = 'inverted_repeats'
        else:
            result['status'] = 'repeats_unclear'

# Write CSV report
with open(args.out_report, 'w', newline='') as csvfile:
    fieldnames = ['contig_id', 'original_length', 'num_hits', 'total_repeat_bp', 
                  'coverage_ratio', 'forward_hits', 'reverse_hits', 'concatemer_detected',
                  'repeat_unit_size', 'num_copies', 'corrected_length', 'status']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(result)

# Write output FASTA (corrected or original)
SeqIO.write([record], args.out_fasta, "fasta")

print(f"Analysis complete. See {args.out_report} for details.")