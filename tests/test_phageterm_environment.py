"""Check runtime command uses the prepared executable, independent of PATH."""
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


class PhageTermEnvironmentTests(unittest.TestCase):
    def test_prepared_executable_and_absolute_inputs_after_chdir(self):
        source=(ROOT/'workflow/rules/3_phageterm.smk').read_text().split('rule phageterm:')[1]
        command=re.search(r'"""(.*?)"""',source,re.S).group(1)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            installation=root/'prepared source'
            binary=installation/'.venv/bin/phageterm'
            binary.parent.mkdir(parents=True)
            (installation/'phagetermvirome').mkdir()
            binary.write_text(f'#!{sys.executable}\nimport json, os, sys\nfrom pathlib import Path\nPath("invocation.json").write_text(json.dumps(dict(argv=sys.argv, cwd=os.getcwd(), pythonpath=os.environ["PYTHONPATH"])))\nPath("Analysis_PhageTerm_report.pdf").write_text("fixture report")\n')
            binary.chmod(0o755)
            trap=root/'wrong_bin'
            trap.mkdir()
            (trap/'phageterm').write_text('#!/bin/sh\nexit 97\n')
            (trap/'phageterm').chmod(0o755)
            (root/'reads.fastq').write_text('fixture reads')
            (root/'virus.fasta').write_text('>a\nACGT\n')
            substitutions={
                '{input.installation:q}':shlex.quote('prepared source'),
                '{input.reads:q}':shlex.quote('reads.fastq'),
                '{input.virus:q}':shlex.quote('virus.fasta'),
                '{output:q}':shlex.quote('sample output/Analysis_PhageTerm_report.pdf'),
                '{log:q}':shlex.quote('run.log'),
                '{threads}':'8',
            }
            for old,new in substitutions.items(): command=command.replace(old,new)
            result=subprocess.run(['bash','-euo','pipefail','-c',command],cwd=root,env=dict(os.environ,PATH=str(trap)+os.pathsep+os.environ['PATH']),capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            invocation=json.loads((root/'sample output/invocation.json').read_text())
            self.assertEqual(invocation['argv'][0],str(binary))
            self.assertEqual(invocation['cwd'],str(root/'sample output'))
            self.assertEqual(invocation['pythonpath'],str(installation/'phagetermvirome'))
            self.assertIn(str(root/'reads.fastq'),invocation['argv'])
            self.assertIn(str(root/'virus.fasta'),invocation['argv'])


if __name__ == '__main__': unittest.main()
