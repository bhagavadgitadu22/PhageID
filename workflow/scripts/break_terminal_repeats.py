#!/usr/bin/env python
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description="Detect and remove second instance of direct terminal repeats from BLAST self-alignment")
    parser.add_argument('--blast', type=str, required=True, help="Path to input BLAST file")
    parser.add_argument('--fasta', type=str, required=True, help="Path to input FASTA file")
    parser.add_argument('--out_report', type=str, required=True, help="Path to output CSV report file")
    parser.add_argument('--out_fasta', type=str, required=True, help="Path to output corrected FASTA file")
    parser.add_argument('--min_identity', type=float, default=80, help="Minimum percent identity (default: 80)")
    parser.add_argument('--max_distance', type=int, default=100, help="Maximum distance from the ends on both sides (default: 20 bp)")
    return parser.parse_args()

args = parse_arguments()

# Parse BLAST results - only one sequence expected
# Collect forward hits above identity threshold
blast_hits = []
with open(args.blast) as f:
    for line in f:
        fields = line.strip().split('\t')
        qseqid, sseqid, pident, length, qstart, qend, qlen, sstart, send, slen = fields
        
        # Skip self-hits at same position
        if qstart == sstart:
            continue
        
        pident = float(pident)
        length = int(length)
        qstart, qend, sstart, send = int(qstart), int(qend), int(sstart), int(send)
        qlen = int(qlen)
        
        # Only keep forward orientation hits above identity threshold
        if pident >= args.min_identity and sstart < send:
            blast_hits.append({
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'length': length,
                'pident': pident,
                'qlen': qlen
            })

# Process the single sequence
record = next(SeqIO.parse(args.fasta, "fasta"))
seq_len = len(record.seq)

# Result dictionary for CSV output
result = {
    'contig_id': record.id,
    'original_length': seq_len,
    'num_hits': 0,
    'dtr_detected': False,
    'dtr_size': 0,
    'dtr_identity': 0.0,
    'five_prime_coords': '',
    'three_prime_coords': '',
    'overlap': 0,
    'cut_position': seq_len,
    'corrected_length': seq_len,
    'status': 'no_hits'
}

if not blast_hits:
    result['status'] = 'no_hits'
else:
    result['num_hits'] = len(blast_hits)
    
    # Find DTR candidates: hits where one copy is at 5' end and other at 3' end
    dtr_candidates = []
    for hit in blast_hits:
        # Check if hit spans both termini within max_distance
        near_5prime = min(hit['qstart'], hit['sstart']) <= args.max_distance
        near_3prime = max(hit['qend'], hit['send']) >= seq_len - args.max_distance
        
        if near_5prime and near_3prime:
            # Determine which is 5' copy and which is 3' copy
            if hit['qstart'] <= args.max_distance:
                # Query is 5' copy, subject is 3' copy
                five_prime_start = hit['qstart']
                five_prime_end = hit['qend']
                three_prime_start = hit['sstart']
                three_prime_end = hit['send']
            else:
                # Subject is 5' copy, query is 3' copy
                five_prime_start = hit['sstart']
                five_prime_end = hit['send']
                three_prime_start = hit['qstart']
                three_prime_end = hit['qend']
            
            dtr_candidates.append({
                'length': hit['length'],
                'pident': hit['pident'],
                'five_prime_start': five_prime_start,
                'five_prime_end': five_prime_end,
                'three_prime_start': three_prime_start,
                'three_prime_end': three_prime_end
            })
    
    if not dtr_candidates:
        result['status'] = 'no_terminal_hits'
    else:
        # Use the longest DTR candidate
        best_dtr = max(dtr_candidates, key=lambda x: x['length'])
        
        result['dtr_detected'] = True
        result['dtr_size'] = best_dtr['length']
        result['dtr_identity'] = best_dtr['pident']
        result['five_prime_coords'] = f"{best_dtr['five_prime_start']}-{best_dtr['five_prime_end']}"
        result['three_prime_coords'] = f"{best_dtr['three_prime_start']}-{best_dtr['three_prime_end']}"
        
        # Calculate cut position: keep entire 5' copy, remove 3' copy
        # If overlap exists, priority is to keep entire 5'
        cut_pos = max(best_dtr['five_prime_end'], best_dtr['three_prime_start'] - 1)
        
        # Check for overlap
        if best_dtr['five_prime_end'] >= best_dtr['three_prime_start']:
            result['overlap'] = best_dtr['five_prime_end'] - best_dtr['three_prime_start'] + 1
        
        result['cut_position'] = cut_pos
        result['corrected_length'] = cut_pos
        result['status'] = 'corrected'
        
        # Truncate sequence to remove 3' DTR
        record.seq = record.seq[:cut_pos]
        record.description = f"{record.description} | DTR removed ({best_dtr['length']} bp, {best_dtr['pident']:.1f}% id)"

# Write CSV report
with open(args.out_report, 'w', newline='') as csvfile:
    fieldnames = ['contig_id', 'original_length', 'num_hits', 'dtr_detected', 
                  'dtr_size', 'dtr_identity', 'five_prime_coords', 'three_prime_coords',
                  'overlap', 'cut_position', 'corrected_length', 'status']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow(result)

# Write output FASTA (corrected or original)
SeqIO.write([record], args.out_fasta, "fasta")

print(f"Analysis complete. See {args.out_report} for details.")