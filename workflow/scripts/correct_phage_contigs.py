#!/usr/bin/env python
import argparse
import os
import subprocess
import tempfile
import shutil
import csv
from pathlib import Path
from multiprocessing import Pool
from Bio import SeqIO


def parse_arguments():
    parser = argparse.ArgumentParser(description="Correct phage contigs by detecting and removing concatemers and DTRs")
    parser.add_argument('--fasta', type=str, required=True, help="Path to input FASTA file with multiple contigs")
    parser.add_argument('--out_fasta', type=str, required=True, help="Path to output corrected FASTA file")
    parser.add_argument('--out_concatemer_report', type=str, required=True, help="Path to output concatemer CSV report")
    parser.add_argument('--out_dtr_report', type=str, required=True, help="Path to output DTR CSV report")
    parser.add_argument('--threads', type=int, default=1, help="Number of threads for parallel processing (default: 1)")
    parser.add_argument('--min_identity', type=float, default=95, help="Minimum percent identity for BLAST (default: 95)")
    parser.add_argument('--min_repeat', type=float, default=90, help="Min repeat threshold for concatemers (default: 90)")
    parser.add_argument('--min_coverage', type=float, default=90, help="Min coverage threshold for concatemers (default: 90)")
    parser.add_argument('--max_distance', type=int, default=20, help="Max distance from ends for DTRs (default: 20)")
    parser.add_argument('--script_dir', type=str, default=None, help="Directory containing break_*.py scripts (default: same as this script)")
    return parser.parse_args()


def run_blast(fasta_file, blast_output):
    """Run BLAST self-alignment on a single contig FASTA file."""
    # Use shell command with proper output redirection
    cmd = (
        f'blastn -query "{fasta_file}" -subject "{fasta_file}" '
        f'-outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen" '
        f'> "{blast_output}"'
    )
    result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"BLAST failed: {result.stderr}")


def process_single_contig(contig_data):
    """Process a single contig through concatemer and DTR correction pipeline."""
    contig_id, seq_record, args, script_dir, temp_base_dir = contig_data
    
    # Create temporary directory for this contig
    temp_dir = tempfile.mkdtemp(dir=temp_base_dir, prefix=f"contig_{contig_id}_")
    
    try:
        # Paths for temporary files
        input_fasta = os.path.join(temp_dir, "input.fasta")
        blast1_output = os.path.join(temp_dir, "blast1.txt")
        concatemer_fasta = os.path.join(temp_dir, "concatemer.fasta")
        concatemer_report = os.path.join(temp_dir, "concatemer.csv")
        blast2_output = os.path.join(temp_dir, "blast2.txt")
        dtr_fasta = os.path.join(temp_dir, "dtr.fasta")
        dtr_report = os.path.join(temp_dir, "dtr.csv")
        
        # Write single contig to file
        SeqIO.write([seq_record], input_fasta, "fasta")
        
        # Step 1: First BLAST
        run_blast(input_fasta, blast1_output)
        
        # Step 2: Break concatemers
        break_concatemers_script = os.path.join(script_dir, "break_concatemers.py")
        cmd_concatemer = [
            'python', break_concatemers_script,
            '--blast', blast1_output,
            '--fasta', input_fasta,
            '--out_report', concatemer_report,
            '--out_fasta', concatemer_fasta,
            '--min_identity', str(args.min_identity),
            '--min_repeat', str(args.min_repeat),
            '--min_coverage', str(args.min_coverage)
        ]
        subprocess.run(cmd_concatemer, check=True, capture_output=True)
        
        # Step 3: Second BLAST on concatemer-corrected sequence
        run_blast(concatemer_fasta, blast2_output)
        
        # Step 4: Break terminal repeats
        break_dtr_script = os.path.join(script_dir, "break_terminal_repeats.py")
        cmd_dtr = [
            'python', break_dtr_script,
            '--blast', blast2_output,
            '--fasta', concatemer_fasta,
            '--out_report', dtr_report,
            '--out_fasta', dtr_fasta,
            '--min_identity', str(args.min_identity),
            '--max_distance', str(args.max_distance)
        ]
        subprocess.run(cmd_dtr, check=True, capture_output=True)
        
        # Read results
        corrected_record = next(SeqIO.parse(dtr_fasta, "fasta"))
        
        with open(concatemer_report, 'r') as f:
            reader = csv.DictReader(f)
            concatemer_row = next(reader)
        
        with open(dtr_report, 'r') as f:
            reader = csv.DictReader(f)
            dtr_row = next(reader)
        
        return {
            'contig_id': contig_id,
            'corrected_record': corrected_record,
            'concatemer_report': concatemer_row,
            'dtr_report': dtr_row,
            'success': True
        }
        
    except Exception as e:
        return {
            'contig_id': contig_id,
            'error': str(e),
            'success': False
        }
    
    finally:
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
        except:
            pass


def main():
    args = parse_arguments()
    
    # Determine script directory
    if args.script_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    else:
        script_dir = args.script_dir
    
    # Load all contigs
    print("Loading contigs from input FASTA...", flush=True)
    contigs = list(SeqIO.parse(args.fasta, "fasta"))
    total_contigs = len(contigs)
    print(f"Loaded {total_contigs} contigs", flush=True)
    
    # Create temporary base directory
    temp_base_dir = tempfile.mkdtemp(prefix="phage_correction_")
    
    try:
        # Prepare data for parallel processing
        contig_data = [
            (record.id, record, args, script_dir, temp_base_dir)
            for record in contigs
        ]
        
        # Process contigs
        print(f"Processing contigs with {args.threads} threads...", flush=True)
        results = []
        
        processed_count = 0
        
        if args.threads > 1:
            with Pool(args.threads) as pool:
                for result in pool.imap_unordered(process_single_contig, contig_data):
                    results.append(result)
                    processed_count += 1
                    if processed_count % 100 == 0 or processed_count == total_contigs:
                        remaining = total_contigs - processed_count
                        print(f"Processed {processed_count}/{total_contigs} contigs. Remaining: {remaining}", flush=True)
        else:
            for data in contig_data:
                result = process_single_contig(data)
                results.append(result)
                processed_count += 1
                if processed_count % 100 == 0 or processed_count == total_contigs:
                    remaining = total_contigs - processed_count
                    print(f"Processed {processed_count}/{total_contigs} contigs. Remaining: {remaining}", flush=True)
        
        # Filter successful results
        successful_results = [r for r in results if r['success']]
        failed_results = [r for r in results if not r['success']]
        
        if failed_results:
            print(f"\nWarning: {len(failed_results)} contigs failed processing:", flush=True)
            for failure in failed_results[:10]:  # Show first 10 failures
                print(f"  - {failure['contig_id']}: {failure['error']}", flush=True)
            if len(failed_results) > 10:
                print(f"  ... and {len(failed_results) - 10} more", flush=True)
        
        print(f"\nSuccessfully processed {len(successful_results)} contigs", flush=True)
        
        # Write combined corrected FASTA
        print("Writing corrected FASTA file...", flush=True)
        corrected_records = [r['corrected_record'] for r in successful_results]
        SeqIO.write(corrected_records, args.out_fasta, "fasta")
        
        # Combine concatemer reports
        print("Writing combined concatemer report...", flush=True)
        concatemer_reports = [r['concatemer_report'] for r in successful_results]
        if concatemer_reports:
            with open(args.out_concatemer_report, 'w', newline='') as f:
                fieldnames = concatemer_reports[0].keys()
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(concatemer_reports)
        
        # Combine DTR reports
        print("Writing combined DTR report...", flush=True)
        dtr_reports = [r['dtr_report'] for r in successful_results]
        if dtr_reports:
            with open(args.out_dtr_report, 'w', newline='') as f:
                fieldnames = dtr_reports[0].keys()
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(dtr_reports)
        
        # Calculate statistics
        print("\n" + "="*60, flush=True)
        print("CORRECTION SUMMARY", flush=True)
        print("="*60, flush=True)
        
        # Concatemer statistics
        concatemers_corrected = sum(1 for r in concatemer_reports if r['status'] == 'corrected')
        total_bp_removed_concatemers = sum(
            int(r['original_length']) - int(r['corrected_length'])
            for r in concatemer_reports if r['status'] == 'corrected'
        )
        print(f"\nConcatemers:", flush=True)
        print(f"  Contigs corrected: {concatemers_corrected}/{len(successful_results)}", flush=True)
        print(f"  Total bp removed: {total_bp_removed_concatemers:,}", flush=True)
        if concatemers_corrected > 0:
            avg_bp = total_bp_removed_concatemers / concatemers_corrected
            print(f"  Average bp removed per corrected contig: {avg_bp:.1f}", flush=True)
        
        # DTR statistics
        dtrs_corrected = sum(1 for r in dtr_reports if r['status'] == 'corrected')
        total_bp_removed_dtrs = sum(
            int(r['original_length']) - int(r['corrected_length'])
            for r in dtr_reports if r['status'] == 'corrected'
        )
        print(f"\nDirect Terminal Repeats:", flush=True)
        print(f"  Contigs corrected: {dtrs_corrected}/{len(successful_results)}", flush=True)
        print(f"  Total bp removed: {total_bp_removed_dtrs:,}", flush=True)
        if dtrs_corrected > 0:
            avg_bp = total_bp_removed_dtrs / dtrs_corrected
            print(f"  Average bp removed per corrected contig: {avg_bp:.1f}", flush=True)
        
        # Overall statistics
        total_corrections = sum(
            1 for r in successful_results 
            if r['concatemer_report']['status'] == 'corrected' or r['dtr_report']['status'] == 'corrected'
        )
        total_bp_removed = total_bp_removed_concatemers + total_bp_removed_dtrs
        print(f"\nOverall:", flush=True)
        print(f"  Total contigs modified: {total_corrections}/{len(successful_results)}", flush=True)
        print(f"  Total bp removed: {total_bp_removed:,}", flush=True)
        
        print("\n" + "="*60, flush=True)
        print(f"Output files:", flush=True)
        print(f"  Corrected FASTA: {args.out_fasta}", flush=True)
        print(f"  Concatemer report: {args.out_concatemer_report}", flush=True)
        print(f"  DTR report: {args.out_dtr_report}", flush=True)
        print("="*60, flush=True)
        
    finally:
        # Clean up temporary base directory
        print("\nCleaning up temporary files...", flush=True)
        try:
            shutil.rmtree(temp_base_dir)
        except:
            pass
    
    print("\nAnalysis complete!", flush=True)


if __name__ == "__main__":
    main()
