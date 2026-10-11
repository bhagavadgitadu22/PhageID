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
    parser.add_argument('--input-rejection', help='Skip assembly after a failed preprocessing attempt')
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
        record.update(Execution_status='completed', Assembly_outcome='assembler_failed')
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
    success = (code == 0 and assembly.is_file() and assembly.stat().st_size > 0
               and all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in args.extra_output))
    # No tool-message parsing: ordinary nonzero exits allow rescue. Signal deaths
    # and missing commands stay incomplete; Slurm status checks catch step OOMs.
    interrupted = code is not None and (code < 0 or 128 < code <= 192)
    execution_error = code is None or code in {126, 127} or interrupted
    record.update(Exit_code=code if code is not None else 'NA',
                  Execution_status='interrupted' if interrupted else 'tool_error' if execution_error else 'completed',
                  Assembly_outcome='not_assessed' if execution_error else 'assembly_produced' if success else 'assembler_failed')
    save()
    if not success:
        assembly.unlink(missing_ok=True)
        for path in args.extra_output:
            Path(path).unlink(missing_ok=True)
        if execution_error:
            raise SystemExit('Assembly execution failed; see ' + str(log_path))
        assembly.write_text('')
        for path in args.extra_output:
            Path(path).write_text('')
    status.write_text('success\n' if success else 'failed\n')


if __name__ == '__main__':
    main()
