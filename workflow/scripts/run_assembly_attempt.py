#!/usr/bin/env python3
"""Record failed assembler commands without blocking the next assembly attempt."""
import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assembly', required=True)
    parser.add_argument('--status', required=True)
    parser.add_argument('--log', required=True)
    parser.add_argument('--extra-output', action='append', default=[])
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('Missing assembler command')
    assembly = Path(args.assembly)
    assembly.parent.mkdir(parents=True, exist_ok=True)
    # Remove an earlier published result before retrying an attempt.
    assembly.unlink(missing_ok=True)
    with open(args.log, 'w') as log:
        try:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            success = (result.returncode == 0 and assembly.is_file() and assembly.stat().st_size > 0
                       and all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in args.extra_output))
            log.write(f'\nAssembly command exit code: {result.returncode}\n')
        except OSError as error:
            success = False
            log.write(f'\nCannot execute assembler: {error}\n')
        if not success:
            log.write('Assembly attempt failed; continuing to rescue or reporting.\n')
    if not success:
        assembly.write_text('')
        for path in args.extra_output:
            Path(path).write_text('')
    Path(args.status).write_text('success\n' if success else 'failed\n')


if __name__ == '__main__':
    main()
