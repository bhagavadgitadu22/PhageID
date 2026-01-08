#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description="Detect and correct concatemers from BLAST self-alignment")
    parser.add_argument('--blast', type=str, required=True, help="Path to BLAST output file")
    parser.add_argument('--fasta', type=str, required=True, help="Path to input FASTA file")
    parser.add_argument('--out_report', type=str, required=True, help="Path to output report file")
    parser.add_argument('--out_fasta', type=str, required=True, help="Path to output corrected FASTA file")
    parser.add_argument('--min_identity', type=float, default=95, help="Minimum percent identity (default: 95)")
    parser.add_argument('--min_repeat', type=float, default=90, help="Repeats are considered when they represent at least x% of the biggest repeat (default: 80)")
    parser.add_argument('--min_coverage', type=float, default=90, help="Minimum coverage percent (default: 80)")
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

with open(args.out_report, 'w') as report:
    report.write(f"Analysis for {record.id} (length: {seq_len} bp)\n")
    report.write(f"{'='*60}\n\n")
    
    if not blast_hits:
        report.write("No significant self-alignments detected.\n")
        report.write("Sequence appears to be a single copy.\n")
    else:
        # Sum total bp in repeats
        total_repeat_bp = sum(hit['length'] for hit in blast_hits)
        coverage_ratio = total_repeat_bp / seq_len
        
        report.write(f"Found {len(blast_hits)} significant self-alignments\n")
        report.write(f"Total repeat bp: {total_repeat_bp} ({coverage_ratio:.1%} of sequence)\n")
        report.write(f"Threshold: {args.min_coverage}%\n\n")
        
        if coverage_ratio < args.min_coverage / 100.0:
            report.write("Coverage below threshold - no concatemer detected.\n")
            report.write("Keeping sequence as-is.\n")
        else:
            report.write("Coverage above threshold - analyzing repeat structure...\n\n")
            
            # Check for inverted repeats (reverse complement)
            reverse_hits = [h for h in blast_hits if h['orientation'] == 'reverse']
            forward_hits = [h for h in blast_hits if h['orientation'] == 'forward']
            
            report.write(f"Forward alignments: {len(forward_hits)}\n")
            report.write(f"Reverse alignments: {len(reverse_hits)}\n\n")
            
            corrected = False
            
            # Check for direct repeats (concatemers)
            if forward_hits:
                # After filtering, all repeats are similar size
                # Use the biggest one as the repeat unit size
                repeat_unit = max(hit['length'] for hit in forward_hits)
                num_copies = round(seq_len / repeat_unit)
                
                if num_copies >= 2:
                    report.write(f"DIRECT REPEAT DETECTED (CONCATEMER)\n")
                    report.write(f"  Repeat unit size: {repeat_unit} bp\n")
                    report.write(f"  Number of copies: ~{num_copies}x\n")
                    report.write(f"  Original length: {seq_len} bp\n")
                    report.write(f"  Corrected length: {repeat_unit} bp\n\n")
                    
                    # Extract single copy
                    record.seq = record.seq[:repeat_unit]
                    record.description = f"{record.description} | corrected from {num_copies}x concatemer"
                    corrected = True
            
            # Check for inverted repeats
            if not corrected and reverse_hits and not forward_hits:
                report.write(f"INVERTED REPEAT DETECTED\n")
                report.write(f"  Sequence contains inverted repeats (reverse complement)\n")
                report.write(f"  This is likely a biological feature, not assembly artifact\n")
                report.write(f"  Keeping sequence as-is\n\n")
            
            if not corrected:
                report.write("Repeat structure unclear or not a simple concatemer.\n")
                report.write("Keeping sequence as-is.\n\n")
                report.write("Alignment details:\n")
                for i, hit in enumerate(blast_hits, 1):
                    report.write(f"  Hit {i}: query {hit['qstart']}-{hit['qend']} vs "
                               f"subject {hit['sstart']}-{hit['send']} "
                               f"({hit['orientation']}, {hit['pident']:.1f}% id, {hit['length']} bp)\n")
    
    # Always write output fasta (corrected or original)
    SeqIO.write([record], args.out_fasta, "fasta")

print(f"Analysis complete. See {args.out_report} for details.")
