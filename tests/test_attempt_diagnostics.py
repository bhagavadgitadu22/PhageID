import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
import assembly_attempt_report as report

class AttemptDiagnosticTests(unittest.TestCase):
    def run_attempt(self, root, program):
        assembly=root/'sample/assembly_attempts/filtered/flye/assembly.fasta'
        log=root/'logs/sample_filtered_flye.log'
        status=assembly.with_name('assembly_status.txt')
        result=subprocess.run([sys.executable,str(ROOT/'workflow/scripts/run_assembly_attempt.py'),
            '--assembly',str(assembly),'--status',str(status),'--log',str(log),'--',
            sys.executable,'-c',program],capture_output=True,
            env=dict(os.environ, SLURM_JOB_ID='123'))
        record=json.loads(Path(str(log)+'.attempt.json').read_text())
        return result, status, record

    def test_killed_command_is_incomplete_and_scheduler_explains_cause(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            result,status,record=self.run_attempt(root,'import os,signal; os.kill(os.getpid(),signal.SIGKILL)')
            self.assertNotEqual(result.returncode,0)
            self.assertFalse(status.exists())
            self.assertEqual(record['Execution_status'],'interrupted')
            self.assertEqual(record['Assembly_outcome'],'not_assessed')
            for state, expected in [('OUT_OF_MEMORY','out_of_memory'),('TIMEOUT','timeout'),('NODE_FAIL','node_failure')]:
                with patch.object(report,'scheduler_states',return_value={'123':state}):
                    report.summarize(root,root/'attempts.tsv')
                with (root/'attempts.tsv').open() as handle:
                    row=next(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual(row['Execution_status'],expected)
                self.assertEqual(row['Sample'],'sample')

    def test_flye_short_read_rejection_allows_rescue(self):
        with tempfile.TemporaryDirectory() as directory:
            result,status,record=self.run_attempt(Path(directory),
                'import sys; print("ERROR: No reads above minimum length threshold (1000)"); sys.exit(1)')
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(status.read_text(),'failed\n')
            self.assertEqual(status.with_name('assembly.fasta').read_text(),'')
            self.assertEqual(record['Execution_status'],'completed')
            self.assertEqual(record['Assembly_outcome'],'assembler_failed')

    def test_flye_no_disjointigs_allows_rescue(self):
        with tempfile.TemporaryDirectory() as directory:
            result,status,record=self.run_attempt(Path(directory),
                'import sys; print("ERROR: No disjointigs were assembled - please check if the read type and genome size parameters are correct"); sys.exit(1)')
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(status.read_text(),'failed\n')
            self.assertEqual(record['Execution_status'],'completed')
            self.assertEqual(record['Assembly_outcome'],'assembler_failed')

    def test_flye_rejection_does_not_mask_interruption(self):
        with tempfile.TemporaryDirectory() as directory:
            result,status,record=self.run_attempt(Path(directory),
                'import sys; print("No reads above minimum length threshold (1000)"); sys.exit(137)')
            self.assertNotEqual(result.returncode,0)
            self.assertFalse(status.exists())
            self.assertEqual(record['Execution_status'],'interrupted')

    def test_ordinary_error_allows_rescue_without_message_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            result,status,record=self.run_attempt(Path(directory),'import sys; sys.exit(1)')
            self.assertEqual(result.returncode,0)
            self.assertEqual(status.read_text(),'failed\n')
            self.assertEqual(record['Execution_status'],'completed')
            self.assertEqual(record['Assembly_outcome'],'assembler_failed')

    def test_zero_exit_without_assembly_allows_rescue(self):
        with tempfile.TemporaryDirectory() as directory:
            result,status,record=self.run_attempt(Path(directory),'pass')
            self.assertEqual(result.returncode,0)
            self.assertEqual(status.read_text(),'failed\n')
            self.assertEqual(record['Assembly_outcome'],'assembler_failed')

    def test_whole_job_interruption_remains_diagnosable(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            _,_,record=self.run_attempt(root,'import sys; sys.exit(1)')
            record['Execution_status']='running'
            (root/'logs/sample_filtered_flye.log.attempt.json').write_text(json.dumps(record))
            with patch.object(report,'scheduler_states',return_value={'123':'TIMEOUT'}):
                report.summarize(root,root/'attempts.tsv')
            with (root/'attempts.tsv').open() as handle:
                row=next(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(row['Execution_status'],'timeout')

    def test_sacct_parser_uses_job_state_not_step_duplicates(self):
        result=subprocess.CompletedProcess([],0,'123|OUT_OF_MEMORY|\n123.batch|FAILED|\n','')
        with patch.object(report.subprocess,'run',return_value=result):
            self.assertEqual(report.scheduler_states({'123'}),{'123':'OUT_OF_MEMORY'})
