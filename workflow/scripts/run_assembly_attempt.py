#!/usr/bin/env python3
"""Keep assessed assembly failures separate from interrupted or broken commands."""
import argparse
import json
import os
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['assembly', 'status', 'log']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--extra-output', action='append', default=[])
    parser.add_argument('--input-rejection', choices=['reads too shallow to subset'])
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command and not args.input_rejection:
        parser.error('Missing assembler command')
    assembly, status, log_path = map(Path, [args.assembly, args.status, args.log])
    assembly.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    assembly.unlink(missing_ok=True)
    status.unlink(missing_ok=True)
    report = Path(str(log_path) + '.attempt.json')
    record = dict(Assembly=str(assembly.resolve()), Log=str(log_path.resolve()),
                  Slurm_job_id=os.environ.get('SLURM_JOB_ID', 'NA'),
                  Execution_status='running', Assembly_outcome='not_assessed', Exit_code='NA')
    def save():
        report.write_text(json.dumps(record, indent=2) + '\n')
    # Persist before execution: whole-job kills cannot run cleanup code.
    save()
    if args.input_rejection:
        log_path.write_text('Assembly input rejected: ' + args.input_rejection + '\n')
        assembly.write_text('')
        status.write_text('failed\n')
        record.update(Execution_status='completed', Assembly_outcome='input_rejected')
        save()
        return
    code = None
    with log_path.open('w') as log:
        try:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            code = result.returncode
            log.write(f'\nAssembly command exit code: {code}\n')
        except OSError as error:
            log.write(f'\nCannot execute assembler: {error}\n')
    text = log_path.read_text(errors='replace')
    success = (code == 0 and assembly.is_file() and assembly.stat().st_size > 0
               and all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in args.extra_output))
    # Exit 20 is reserved by our consensus script for completed QC with no usable cluster.
    # Known SPAdes/Flye read rejections are assessed input limitations.
    qc_failed = code == 20 and 'No clusters passed Autocycler QC.' in text
    input_rejected = ((code == 64 and 'Invalid kmer coverage histogram' in text)
                      or (code == 1 and 'No reads above minimum length threshold (' in text)
                      or (code == 21 and 'Autocycler requires at least two successful input assemblies.' in text))
    if 'SIGKILL' in text:
        input_rejected = qc_failed = False
    assessed_failure = qc_failed or input_rejected
    record.update(Exit_code=code if code is not None else 'NA',
                  Execution_status='completed' if success or assessed_failure else 'interrupted' if code == -9 or 'SIGKILL' in text else 'tool_error',
                  Assembly_outcome='assembly_produced' if success else 'qc_failed' if qc_failed
                  else 'input_rejected' if input_rejected else 'not_assessed')
    save()
    if not success:
        assembly.unlink(missing_ok=True)
        for path in args.extra_output:
            Path(path).unlink(missing_ok=True)
        if not assessed_failure:
            raise SystemExit('Assembly execution failed; see ' + str(log_path))
        assembly.write_text('')
        for path in args.extra_output:
            Path(path).write_text('')
    status.write_text('success\n' if success else 'failed\n')


if __name__ == '__main__':
    main()
