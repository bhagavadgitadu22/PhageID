import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / 'workflow/scripts/run_assembly_attempt.py'


class AssemblyFailureTests(unittest.TestCase):
    def test_failed_command_discards_partial_assembly_and_keeps_log(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            assembly, status, log, info = (p / name for name in ['assembly.fa', 'status', 'log', 'info'])
            assembly.write_text('stale result')
            command = [sys.executable, str(WRAPPER), '--assembly', str(assembly),
                       '--status', str(status), '--log', str(log), '--extra-output', str(info), '--',
                       sys.executable, '-c',
                       'import pathlib,sys; pathlib.Path(sys.argv[1]).write_text(">partial\\nACGT\\n"); print("coverage error"); sys.exit(64)', str(assembly)]
            subprocess.run(command, check=True)
            self.assertEqual(assembly.read_text(), '')
            self.assertEqual(info.read_text(), '')
            self.assertEqual(status.read_text(), 'failed\n')
            self.assertIn('coverage error', log.read_text())
            self.assertIn('exit code: 64', log.read_text())
            command[-2] = 'import pathlib,sys; pathlib.Path(sys.argv[1]).write_text(">c\\nACGT\\n"); pathlib.Path(sys.argv[1]).with_name("info").write_text("metadata")'
            subprocess.run(command, check=True)
            self.assertEqual(status.read_text(), 'success\n')
            self.assertEqual(assembly.read_text(), '>c\nACGT\n')

    def test_missing_executable_is_failed_attempt(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            subprocess.run([sys.executable, str(WRAPPER), '--assembly', str(p/'assembly'),
                            '--status', str(p/'status'), '--log', str(p/'log'), '--',
                            str(p/'missing_assembler')], check=True)
            self.assertEqual((p/'status').read_text(), 'failed\n')
            self.assertIn('Cannot execute assembler', (p/'log').read_text())


if __name__ == '__main__':
    unittest.main()
