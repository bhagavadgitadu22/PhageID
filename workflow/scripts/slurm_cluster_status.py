#!/usr/bin/env python3
"""Snakemake cluster-status command: resource failures in any job step win."""
import argparse
import re
import subprocess

FAILURES = {'OUT_OF_MEMORY', 'TIMEOUT', 'NODE_FAIL', 'BOOT_FAIL', 'CANCELLED',
            'PREEMPTED', 'FAILED', 'DEADLINE', 'REVOKED'}
PRIORITY = ['OUT_OF_MEMORY', 'TIMEOUT', 'NODE_FAIL', 'BOOT_FAIL', 'CANCELLED',
            'PREEMPTED', 'DEADLINE', 'REVOKED', 'FAILED']

def accounting_states(ids):
    result = subprocess.run(['sacct', '--jobs', ','.join(sorted(ids)), '--noheader',
                             '--parsable2', '--format=JobIDRaw,State'],
                            capture_output=True, text=True, timeout=30, check=True)
    primary, steps = {}, {}
    for line in result.stdout.splitlines():
        fields = line.split('|')
        if len(fields) < 2 or not fields[1].strip():
            continue
        job, state = fields[0].strip(), fields[1].split()[0].rstrip('+')
        parent = job.split('.')[0]
        if parent not in ids:
            continue
        steps.setdefault(parent, set()).add(state)
        if job == parent:
            primary[parent] = state
    for job, states in steps.items():
        failure = next((state for state in PRIORITY if state in states), None)
        if failure:
            primary[job] = failure
    return primary

def job_status(job_id):
    state = accounting_states({job_id}).get(job_id)
    if state in FAILURES:
        return 'failed'
    if state == 'COMPLETED':
        return 'success'
    # Accounting can lag behind submission; unknown state is not success.
    return 'running'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job_id')
    args = parser.parse_args()
    # Existing sbatch templates return "Submitted batch job N" instead of --parsable.
    match = re.fullmatch(r'(?:Submitted batch job\s+)?(\d+)', args.job_id.strip())
    if not match:
        parser.error('Expected a Slurm job ID or "Submitted batch job N"')
    print(job_status(match.group(1)))
