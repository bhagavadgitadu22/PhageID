import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1]
SOURCE=(ROOT/'workflow/rules/02_assembly.smk').read_text()


class AutocyclerTests(unittest.TestCase):
    def test_consensus_steps_and_qc_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            executable=p/'autocycler'
            executable.write_text(r"""#!/usr/bin/env python3
import os,sys
from pathlib import Path
args=sys.argv[1:]
step=args[0]
with open(os.environ['CALLS'],'a') as log: log.write(step+'\n')
if step=='compress':
    Path(args[args.index('-a')+1]).mkdir()
if step=='cluster' and os.environ.get('FAIL_QC')!='1':
    (Path(args[args.index('-a')+1])/'clustering/qc_pass/cluster_001').mkdir(parents=True)
if step=='resolve':
    (Path(args[args.index('-c')+1])/'5_final.gfa').write_text('graph')
if step=='combine':
    (Path(args[args.index('-a')+1])/'consensus_assembly.fasta').write_text('>consensus\nACGT\n')
""")
            executable.chmod(0o755)
            command=SOURCE.split('rule assembly_reads_autocycler:')[1].split('"""')[1]
            quote=lambda value: shlex.quote(str(value))
            command=command.replace(':q}', '}').format(
                params=SimpleNamespace(work=quote(p/'work')),
                input=SimpleNamespace(reads=quote(p/'reads.fastq'),candidates=quote(p/'candidates')),
                output=quote(p/'assembly.fa'),log=quote(p/'run.log'),threads=2)
            env=dict(os.environ,PATH=str(p)+os.pathsep+os.environ['PATH'],CALLS=str(p/'calls'))
            subprocess.run(['bash','-euo','pipefail','-c',command],env=env,check=True)
            self.assertEqual((p/'calls').read_text().splitlines(),['compress','cluster','trim','resolve','combine'])
            self.assertEqual((p/'assembly.fa').read_text(),'>consensus\nACGT\n')
            (p/'assembly.fa').unlink()
            env['FAIL_QC']='1'
            result=subprocess.run(['bash','-euo','pipefail','-c',command],env=env)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((p/'assembly.fa').exists())

    def test_candidate_failures_keep_successful_assemblies(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            subsets=p/'subsets';subsets.mkdir()
            (subsets/'sample_01.fastq').write_text('@r\nACGT\n+\nIIII\n')
            (p/'size').write_text('50000\n')
            executable=p/'autocycler'
            executable.write_text(r"""#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
if args[1]=='flye': sys.exit(1)
Path(args[args.index('--out_prefix')+1]+'.fasta').write_text('>contig\nACGT\n')
""")
            executable.chmod(0o755)
            command=SOURCE.split('rule autocycler_candidate_assemblies:')[1].split('"""')[1]
            quote=lambda value: shlex.quote(str(value))
            command=command.replace(':q}', '}').format(
                input=SimpleNamespace(subsets=quote(subsets),genome_size=quote(p/'size')),
                output=quote(p/'candidates'),params=SimpleNamespace(read_type='ont_r10'),threads=2,log=quote(p/'run.log'))
            env=dict(os.environ,PATH=str(p)+os.pathsep+os.environ['PATH'])
            subprocess.run(['bash','-euo','pipefail','-c',command],env=env,check=True)
            self.assertEqual(len(list((p/'candidates/assemblies').glob('*.fasta'))),7)
            self.assertIn('Assembly failed: flye',(p/'run.log').read_text())
