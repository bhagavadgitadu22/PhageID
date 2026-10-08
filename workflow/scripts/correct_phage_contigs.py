#!/usr/bin/env python3
"""Correct likely viral concatemers and direct terminal repeats."""

import argparse
import copy
import csv
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

from Bio import SeqIO

CONCATEMER_FIELDS = (
    "contig_id", "original_length", "num_hits", "unique_repeat_coverage_bp",
    "unique_repeat_coverage_ratio", "forward_hits", "reverse_hits", "concatemer_detected",
    "repeat_unit_size", "num_copies", "corrected_length", "status",
)
DTR_FIELDS = (
    "contig_id", "original_length", "num_hits", "dtr_detected", "dtr_size",
    "dtr_identity", "five_prime_coords", "three_prime_coords", "overlap",
    "cut_position", "corrected_length", "status",
)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="Input multi-contig FASTA")
    parser.add_argument("--out_fasta", required=True, help="Corrected FASTA")
    parser.add_argument("--out_concatemer_report", required=True, help="Concatemer CSV")
    parser.add_argument("--out_dtr_report", required=True, help="DTR CSV")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min_identity", type=float, default=95)
    parser.add_argument("--min_repeat", type=float, default=90)
    parser.add_argument("--min_coverage", type=float, default=90)
    parser.add_argument("--max_distance", type=int, default=20)
    return parser.parse_args()

def run_blast(fasta_file, blast_output):
    """Run BLAST self-alignment on one contig FASTA."""
    command = [
        "blastn", "-query", str(fasta_file), "-subject", str(fasta_file),
        "-outfmt", "6 qseqid sseqid pident length qstart qend qlen sstart send slen",
    ]
    with open(blast_output, "w") as output_handle:
        subprocess.run(
            command, check=True, stdout=output_handle,
            stderr=subprocess.PIPE, text=True,
        )

def read_blast_hits(path):
    """Read the fixed ten-column self-BLAST format once."""
    hits = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ValueError(f"Expected 10 BLAST columns in {path}:{line_number}")
            qseqid, sseqid, pident, length, qstart, qend, qlen, sstart, send, slen = fields
            hits.append({
                "qseqid": qseqid, "sseqid": sseqid, "pident": float(pident),
                "length": int(length), "qstart": int(qstart), "qend": int(qend),
                "qlen": int(qlen), "sstart": int(sstart), "send": int(send),
                "slen": int(slen),
            })
    return hits

def merged_query_coverage_bp(hits):
    """Return unique bp covered by inclusive BLAST query intervals."""
    intervals = sorted(
        (min(hit["qstart"], hit["qend"]), max(hit["qstart"], hit["qend"]))
        for hit in hits
    )
    if not intervals:
        return 0
    covered = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    return covered + current_end - current_start + 1

def correct_concatemer(record, hits, min_identity, min_repeat, min_coverage):
    """Infer a repeat period from self-alignment offsets and verify every copy."""
    sequence_length = len(record.seq)
    result = {
        "contig_id": record.id, "original_length": sequence_length,
        "num_hits": 0, "unique_repeat_coverage_bp": 0,
        "unique_repeat_coverage_ratio": 0.0,
        "forward_hits": 0, "reverse_hits": 0, "concatemer_detected": False,
        "repeat_unit_size": 0, "num_copies": 0,
        "corrected_length": sequence_length, "status": "no_repeats",
    }
    eligible = [
        {**hit, "orientation": "forward" if hit["sstart"] < hit["send"] else "reverse"}
        for hit in hits
        if hit["qstart"] != hit["sstart"] and hit["pident"] >= min_identity
    ]
    if eligible:
        largest = max(hit["length"] for hit in eligible)
        eligible = [hit for hit in eligible if hit["length"] >= largest * min_repeat / 100.0]
    if not eligible:
        return record, result

    forward = [hit for hit in eligible if hit["orientation"] == "forward"]
    reverse = [hit for hit in eligible if hit["orientation"] == "reverse"]
    unique_coverage_bp = merged_query_coverage_bp(eligible + [
        {"qstart": hit["sstart"], "qend": hit["send"]} for hit in eligible
    ])
    coverage_ratio = unique_coverage_bp / sequence_length
    result.update({
        "num_hits": len(eligible), "unique_repeat_coverage_bp": unique_coverage_bp,
        "unique_repeat_coverage_ratio": coverage_ratio, "forward_hits": len(forward),
        "reverse_hits": len(reverse),
    })
    if coverage_ratio < min_coverage / 100.0:
        result["status"] = "below_threshold"
    elif forward:
        result["status"] = "repeats_unclear"
        sequence = str(record.seq).upper()
        # A 1..200 vs 101..300 alignment supports a 100-bp period,
        # even though the alignment itself is 200 bp long.
        offsets = sorted({abs(hit["sstart"] - hit["qstart"]) for hit in forward
                          if hit["qstart"] < hit["qend"]
                          and hit["sstart"] - hit["qstart"] == hit["send"] - hit["qend"]})
        for repeat_unit in offsets:
            copies, remainder = divmod(sequence_length, repeat_unit)
            # Be conservative: partial copies or indels need more evidence.
            if copies < 2 or remainder:
                continue
            supporting_hits = [hit for hit in forward
                               if abs(hit["sstart"] - hit["qstart"]) == repeat_unit
                               and hit["sstart"] - hit["qstart"] == hit["send"] - hit["qend"]]
            # Include both sides of each alignment; reciprocal hits are optional.
            intervals = supporting_hits + [
                {"qstart": hit["sstart"], "qend": hit["send"]} for hit in supporting_hits
            ]
            if merged_query_coverage_bp(intervals) / sequence_length < min_coverage / 100.0:
                continue
            unit = sequence[:repeat_unit]
            if any(
                sum(a == b and a in "ACGT" for a, b in zip(unit, sequence[start:start + repeat_unit]))
                / repeat_unit < min_identity / 100.0
                for start in range(repeat_unit, sequence_length, repeat_unit)
            ):
                continue
            result.update({
                "concatemer_detected": True, "repeat_unit_size": repeat_unit,
                "num_copies": copies, "corrected_length": repeat_unit,
                "status": "corrected",
            })
            record.seq = record.seq[:repeat_unit]
            record.description = f"{record.description} | corrected from {copies}x concatemer"
            break
    elif reverse:
        result["status"] = "inverted_repeats"
    else:
        result["status"] = "repeats_unclear"
    return record, result

def correct_dtr(record, hits, min_identity, max_distance):
    """Remove the second copy of the longest qualifying direct terminal repeat."""
    sequence_length = len(record.seq)
    result = {
        "contig_id": record.id, "original_length": sequence_length,
        "num_hits": 0, "dtr_detected": False, "dtr_size": 0,
        "dtr_identity": 0.0, "five_prime_coords": "",
        "three_prime_coords": "", "overlap": 0,
        "cut_position": sequence_length, "corrected_length": sequence_length,
        "status": "no_hits",
    }
    eligible = [
        hit for hit in hits
        if hit["qstart"] != hit["sstart"]
        and hit["pident"] >= min_identity
        and hit["sstart"] < hit["send"]
    ]
    if not eligible:
        return record, result
    result["num_hits"] = len(eligible)
    candidates = []
    for hit in eligible:
        near_five_prime = min(hit["qstart"], hit["sstart"]) <= max_distance
        near_three_prime = max(hit["qend"], hit["send"]) >= sequence_length - max_distance
        if not (near_five_prime and near_three_prime):
            continue
        if hit["qstart"] <= max_distance:
            five_start, five_end = hit["qstart"], hit["qend"]
            three_start, three_end = hit["sstart"], hit["send"]
        else:
            five_start, five_end = hit["sstart"], hit["send"]
            three_start, three_end = hit["qstart"], hit["qend"]
        candidates.append({
            "length": hit["length"], "pident": hit["pident"],
            "five_prime_start": five_start, "five_prime_end": five_end,
            "three_prime_start": three_start, "three_prime_end": three_end,
        })
    if not candidates:
        result["status"] = "no_terminal_hits"
        return record, result

    best = max(candidates, key=lambda hit: hit["length"])
    cut_position = max(best["five_prime_end"], best["three_prime_start"] - 1)
    overlap = max(0, best["five_prime_end"] - best["three_prime_start"] + 1)
    result.update({
        "dtr_detected": True, "dtr_size": best["length"],
        "dtr_identity": best["pident"],
        "five_prime_coords": f"{best['five_prime_start']}-{best['five_prime_end']}",
        "three_prime_coords": f"{best['three_prime_start']}-{best['three_prime_end']}",
        "overlap": overlap, "cut_position": cut_position,
        "corrected_length": cut_position, "status": "corrected",
    })
    record.seq = record.seq[:cut_position]
    record.description = (
        f"{record.description} | DTR removed "
        f"({best['length']} bp, {best['pident']:.1f}% id)"
    )
    return record, result

def exception_text(error):
    message = str(error)
    if isinstance(error, subprocess.CalledProcessError) and error.stderr:
        stderr = error.stderr.decode() if isinstance(error.stderr, bytes) else error.stderr
        message = f"{message}; stderr: {stderr.strip()}"
    return message

def staged_output(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", prefix=f".{destination.name}.", suffix=".tmp",
        dir=destination.parent, delete=False,
    )
    handle.close()
    return Path(handle.name)

def publish_empty_outputs(args):
    destinations = (
        Path(args.out_fasta), Path(args.out_concatemer_report), Path(args.out_dtr_report),
    )
    staged = [staged_output(destination) for destination in destinations]
    try:
        for path, fields in zip(staged[1:], (CONCATEMER_FIELDS, DTR_FIELDS)):
            with path.open("w", newline="") as handle:
                csv.DictWriter(handle, fieldnames=fields).writeheader()
        for source, destination in zip(staged, destinations):
            os.replace(source, destination)
        staged = []
    finally:
        for path in staged:
            path.unlink(missing_ok=True)

def validate_result(input_record, result):
    contig_id = input_record.id
    corrected = result["corrected_record"]
    concatemer = result["concatemer_report"]
    dtr = result["dtr_report"]
    if corrected.id != contig_id:
        raise ValueError(f"Corrected FASTA ID changed from {contig_id!r} to {corrected.id!r}")
    if concatemer.get("contig_id") != contig_id or dtr.get("contig_id") != contig_id:
        raise ValueError(f"Correction report ID mismatch for {contig_id!r}")
    if not corrected.seq:
        raise ValueError(f"Correction produced an empty sequence for {contig_id!r}")
    original_length = len(input_record.seq)
    transitions = (
        ("concatemer", int(concatemer["original_length"]), int(concatemer["corrected_length"])),
        ("DTR", int(dtr["original_length"]), int(dtr["corrected_length"])),
    )
    if transitions[0][1] != original_length:
        raise ValueError(f"Concatemer input length mismatch for {contig_id!r}")
    if transitions[1][1] != transitions[0][2]:
        raise ValueError(f"DTR input length mismatch for {contig_id!r}")
    if transitions[1][2] != len(corrected.seq):
        raise ValueError(f"Final corrected length mismatch for {contig_id!r}")
    for (stage, before, after), row in zip(transitions, (concatemer, dtr)):
        if after != before and row.get("status") != "corrected":
            raise ValueError(f"{stage} length changed without a corrected report for {contig_id!r}")

def process_single_contig(contig_data):
    contig_id, record, args, temp_base_dir = contig_data
    original_length = len(record.seq)
    record = copy.deepcopy(record)
    temp_dir = tempfile.mkdtemp(dir=temp_base_dir, prefix="contig_")
    stage = "initialization"
    try:
        input_fasta = Path(temp_dir) / "input.fasta"
        first_blast = Path(temp_dir) / "blast1.tsv"
        second_blast = Path(temp_dir) / "blast2.tsv"
        stage = "writing input FASTA"
        SeqIO.write([record], input_fasta, "fasta")
        stage = "initial self-BLAST"
        run_blast(input_fasta, first_blast)
        stage = "concatemer correction"
        record, concatemer = correct_concatemer(record, read_blast_hits(first_blast), args.min_identity,args.min_repeat, args.min_coverage)
        stage = "writing concatemer-corrected FASTA"
        SeqIO.write([record], input_fasta, "fasta")
        stage = "post-concatemer self-BLAST"
        run_blast(input_fasta, second_blast)
        stage = "DTR correction"
        record, dtr = correct_dtr(record, read_blast_hits(second_blast), args.min_identity, args.max_distance)
        return {
            "contig_id": contig_id, "corrected_record": record,
            "concatemer_report": concatemer, "dtr_report": dtr, "success": True,
        }
    except Exception as error:
        return {
            "contig_id": contig_id, "stage": stage,
            "original_length": original_length,
            "error": exception_text(error), "success": False,
        }
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def write_outputs(args, input_records, results):
    destinations = (args.out_fasta, args.out_concatemer_report, args.out_dtr_report)
    staged = [staged_output(destination) for destination in destinations]
    try:
        if SeqIO.write([result["corrected_record"] for result in results], staged[0], "fasta") != len(results):
            raise ValueError("Corrected FASTA record count mismatch")
        for path, fields, key in (
            (staged[1], CONCATEMER_FIELDS, "concatemer_report"),
            (staged[2], DTR_FIELDS, "dtr_report"),
        ):
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(result[key] for result in results)
        for input_record, result in zip(input_records, results):
            validate_result(input_record, result)
        for path in staged:
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Staged output is missing or empty: {path}")
        for source, destination in zip(staged, destinations):
            os.replace(source, destination)
        staged = []
    finally:
        for path in staged:
            path.unlink(missing_ok=True)

def main():
    args = parse_arguments()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    print("Loading contigs from input FASTA...", flush=True)
    contigs = list(SeqIO.parse(args.fasta, "fasta"))
    if not contigs:
        print(f"Input FASTA contains no viral sequences: {args.fasta}", flush=True)
        publish_empty_outputs(args)
        return
    input_ids = [record.id for record in contigs]
    duplicates = sorted(name for name, count in Counter(input_ids).items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate input contig IDs: {duplicates[:10]}")
    print(f"Loaded {len(contigs)} contigs", flush=True)

    temp_base_dir = tempfile.mkdtemp(prefix="phage_correction_")
    try:
        work = [(record.id, record, args, temp_base_dir) for record in contigs]
        if args.threads > 1:
            with Pool(args.threads) as pool:
                results = list(pool.imap_unordered(process_single_contig, work))
        else:
            results = [process_single_contig(item) for item in work]
    finally:
        shutil.rmtree(temp_base_dir, ignore_errors=True)

    failures = [result for result in results if not result["success"]]
    if failures:
        for failure in sorted(failures, key=lambda item: item["contig_id"]):
            print(
                f"ERROR: {failure['contig_id']} (stage={failure['stage']}): "
                f"{failure['error']}", flush=True,
            )
        raise RuntimeError(
            f"Correction failed for {len(failures)}/{len(contigs)} contigs; no outputs published"
        )
    by_id = {result["contig_id"]: result for result in results}
    if set(by_id) != set(input_ids) or len(results) != len(contigs):
        raise ValueError("Correction result ID/count mismatch")
    ordered = [by_id[contig_id] for contig_id in input_ids]
    write_outputs(args, contigs, ordered)

    concatemer_count = sum(r["concatemer_report"]["status"] == "corrected" for r in ordered)
    dtr_count = sum(r["dtr_report"]["status"] == "corrected" for r in ordered)
    concatemer_bp = sum(
        r["concatemer_report"]["original_length"] - r["concatemer_report"]["corrected_length"]
        for r in ordered if r["concatemer_report"]["status"] == "corrected"
    )
    dtr_bp = sum(
        r["dtr_report"]["original_length"] - r["dtr_report"]["corrected_length"]
        for r in ordered if r["dtr_report"]["status"] == "corrected"
    )
    print(
        f"Corrected {concatemer_count} concatemers ({concatemer_bp} bp removed) and "
        f"{dtr_count} DTRs ({dtr_bp} bp removed) across {len(contigs)} contigs.",
        flush=True,
    )

if __name__ == "__main__":
    main()
