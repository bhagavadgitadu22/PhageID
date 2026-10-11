import os
import tempfile
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
import slurm_cluster_status as status


class SlurmStatusTests(unittest.TestCase):
    def check(self, output):
        result=subprocess.CompletedProcess([],0,output,'')
        with patch.object(status.subprocess,'run',return_value=result):
            return status.job_status('123')

    def test_completed_job_with_oom_step_fails(self):
        self.assertEqual(self.check('123|COMPLETED|\n123.batch|OUT_OF_MEMORY|\n'),'failed')

    def test_timeout_and_node_failures_fail(self):
        for state in ['TIMEOUT','NODE_FAIL','CANCELLED by 456','PREEMPTED']:
            self.assertEqual(self.check('123|'+state+'|\n'),'failed')

    def test_completed_succeeds_and_missing_accounting_waits(self):
        self.assertEqual(self.check('123|COMPLETED|\n123.batch|COMPLETED|\n'),'success')
        self.assertEqual(self.check(''),'running')
        self.assertEqual(self.check('123|RUNNING|\n'),'running')

    def test_cli_accepts_sbatch_output_and_flags_oom_step(self):
        with tempfile.TemporaryDirectory() as directory:
            executable=Path(directory)/'sacct'
            executable.write_text('#!/bin/sh\nprintf "123|COMPLETED|\\n123.batch|OUT_OF_MEMORY|\\n"\n')
            executable.chmod(0o755)
            result=subprocess.run([sys.executable,status.__file__,'Submitted batch job 123'],
                capture_output=True,text=True,check=True,
                env=dict(os.environ,PATH=directory+os.pathsep+os.environ['PATH']))
            self.assertEqual(result.stdout.strip(),'failed')

    def test_accounting_errors_do_not_report_success(self):
        with patch.object(status.subprocess,'run',side_effect=subprocess.CalledProcessError(1,['sacct'])):
            with self.assertRaises(subprocess.CalledProcessError):
                status.job_status('123')
