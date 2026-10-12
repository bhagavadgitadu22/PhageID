#!/usr/bin/env python3
"""Summarize assembly attempts and query Slurm for interrupted jobs."""
import argparse
import csv
import json
import subprocess
from pathlib import Path
from slurm_cluster_status import accounting_states

FIELDS = ['Sample', 'Readset', 'Assembler', 'Subset', 'Execution_status', 'Assembly_outcome',
          'Exit_code', 'Slurm_job_id', 'Slurm_state', 'Log']
FAILURES = {'OUT_OF_MEMORY': 'out_of_memory', 'TIMEOUT': 'timeout', 'NODE_FAIL': 'node_failure',
            'CANCELLED': 'cancelled', 'PREEMPTED': 'preempted', 'BOOT_FAIL': 'node_failure',
            'FAILED': 'execution_failed', 'DEADLINE': 'timeout', 'REVOKED': 'cancelled'}

def scheduler_states(ids):
    if not ids:
        return {}
    try:
        return accounting_states(ids)
    except (OSError, subprocess.SubprocessError):
        return {}

def summarize(results, output):
    results = Path(results).resolve()
    records = [json.loads(path.read_text()) for path in sorted((results/'logs').glob('*.attempt.json'))]
    states = scheduler_states({str(row['Slurm_job_id']) for row in records if row['Slurm_job_id'] != 'NA'})
    rows, recorded = [], set()
    for record in records:
        assembly = Path(record['Assembly'])
        parts = assembly.relative_to(results).parts
        if len(parts) < 5 or parts[1] != 'assembly_attempts':
            raise ValueError(f'Unexpected attempt path: {assembly}')
        sample, _, readset, assembler = parts[:4]
        subset = parts[-2] if 'candidate_jobs' in parts else 'NA'
        if subset != 'NA':
            assembler = parts[-3]
        recorded.add(assembly)
        state = states.get(str(record['Slurm_job_id']), 'unavailable')
        execution, outcome = record['Execution_status'], record['Assembly_outcome']
        if state in FAILURES:
            execution, outcome = FAILURES[state], 'not_assessed'
        elif execution == 'running' and state in {'COMPLETED', 'unavailable'}:
            execution = 'assessment_incomplete'
        rows.append(dict(zip(FIELDS, [sample, readset, assembler, subset, execution, outcome,
                                      record['Exit_code'], record['Slurm_job_id'], state, record['Log']])))
    for status in sorted(results.glob('*/assembly_attempts/*/*/assembly_status.txt')):
        if (status.parent/'assembly.fasta').resolve() in recorded:
            continue
        sample, _, readset, assembler, _ = status.relative_to(results).parts
        success = status.read_text().strip() == 'success'
        rows.append(dict(zip(FIELDS, [sample, readset, assembler, 'NA', 'legacy_unknown',
                                      'assembly_produced' if success else 'not_assessed', 'NA', 'NA', 'unavailable', 'NA'])))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    summarize(args.results_dir, args.output)
