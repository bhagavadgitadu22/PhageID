import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1]
SOURCE='\n'.join(path.read_text() for path in sorted((ROOT/'workflow/rules').glob('*.smk')))


class AutocyclerTests(unittest.TestCase):
    def test_shallow_subsampling_skips_candidates_and_allows_rescue(self):
        import json
        import sys
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/'size').write_text('100000\n')
            executable=p/'autocycler'
            executable.write_text('#!/bin/sh\necho "Error: input reads are too shallow to subset"\nexit 1\n')
            executable.chmod(0o755)
            env=dict(os.environ,PATH=str(p)+os.pathsep+os.environ['PATH'])
            quote=lambda value: shlex.quote(str(value))
            command=SOURCE.split('rule autocycler_subsample:')[1].split('"""')[1]
            command=command.replace(':q}', '}').format(
                input=SimpleNamespace(reads=quote(p/'reads.fastq'),genome_size=quote(p/'size')),
                output=quote(p/'subsets'),params=SimpleNamespace(subset_count=4),log=quote(p/'subsample.log'))
            subprocess.run(['bash','-euo','pipefail','-c',command],env=env,check=True)
            self.assertTrue((p/'subsets/input_rejected.txt').exists())
            candidate=SOURCE.split('rule autocycler_candidate_assembly:')[1].split('"""')[1]
            candidate=candidate.replace('./scripts/',str(ROOT/'workflow/scripts')+'/')
            candidate=candidate.replace(':q}', '}').format(
                params=SimpleNamespace(plassembler_db=quote(p/'db'),read_type='ont_r10',reads=quote(p/'subsets/sample_01.fastq')),
                input=SimpleNamespace(genome_size=quote(p/'size')),output=quote(p/'candidate'),
                wildcards=SimpleNamespace(candidate_assembler='flye'),threads=1,log=quote(p/'candidate.log'))
            subprocess.run(['bash','-euo','pipefail','-c',candidate],env=env,check=True)
            self.assertEqual((p/'candidate/status.txt').read_text(),'failed\n')
            self.assertEqual(json.loads((p/'candidate.log.attempt.json').read_text())['Assembly_outcome'],'assembler_failed')
            # No candidates is an assessed input limitation, so selection can try all reads.
            (p/'candidates/assemblies').mkdir(parents=True)
            subprocess.run([sys.executable,str(ROOT/'workflow/scripts/run_assembly_attempt.py'),
                '--assembly',str(p/'consensus.fasta'),'--status',str(p/'status'),'--log',str(p/'consensus.log'),'--',
                'bash',str(ROOT/'workflow/scripts/autocycler_consensus.sh'),str(p/'work'),str(p/'candidates'),str(p/'reads.fastq'),'1',str(p/'consensus.fasta')],check=True)
            self.assertEqual((p/'status').read_text(),'failed\n')
            self.assertEqual(json.loads((p/'consensus.log.attempt.json').read_text())['Assembly_outcome'],'assembler_failed')
            # Ordinary subsampling errors also allow rescue without message parsing.
            executable.write_text('#!/bin/sh\necho "Unexpected error"\nexit 1\n')
            result=subprocess.run(['bash','-euo','pipefail','-c',command],env=env)
            self.assertEqual(result.returncode,0)
            self.assertTrue((p/'subsets/input_rejected.txt').exists())
            executable.write_text('#!/bin/sh\nexit 137\n')
            result=subprocess.run(['bash','-euo','pipefail','-c',command],env=env)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((p/'subsets/input_rejected.txt').exists())

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
            (p/'candidates/assemblies').mkdir(parents=True)
            for name in ['a','b']: (p/'candidates/assemblies'/f'{name}.fasta').write_text('>c\nACGT\n')
            command=SOURCE.split('rule assembly_reads_autocycler:')[1].split('"""')[1]
            command=command.replace('./scripts/',str(ROOT/'workflow/scripts')+'/')
            quote=lambda value: shlex.quote(str(value))
            command=command.replace(':q}', '}').format(
                params=SimpleNamespace(work=quote(p/'work')),
                input=SimpleNamespace(reads=quote(p/'reads.fastq'),candidates=quote(p/'candidates')),
                output=SimpleNamespace(assembly=quote(p/'assembly.fa'),status=quote(p/'status')),log=quote(p/'run.log'),threads=2)
            env=dict(os.environ,PATH=str(p)+os.pathsep+os.environ['PATH'],CALLS=str(p/'calls'))
            subprocess.run(['bash','-euo','pipefail','-c',command],env=env,check=True)
            self.assertEqual((p/'calls').read_text().splitlines(),['compress','cluster','trim','resolve','combine'])
            self.assertEqual((p/'assembly.fa').read_text(),'>consensus\nACGT\n')
            (p/'assembly.fa').unlink()
            env['FAIL_QC']='1'
            result=subprocess.run(['bash','-euo','pipefail','-c',command],env=env)
            self.assertEqual(result.returncode,0)
            self.assertEqual((p/'status').read_text().strip(),'failed')
            self.assertEqual((p/'assembly.fa').read_text(),'')

    def test_candidate_ordinary_failures_allow_consensus(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            subsets=p/'subsets';subsets.mkdir()
            (subsets/'sample_01.fastq').write_text('@r\nACGT\n+\nIIII\n')
            (p/'size').write_text('50000\n')
            executable=p/'autocycler'
            executable.write_text(r"""#!/usr/bin/env python3
import os,sys
from pathlib import Path
args=sys.argv[1:]
if args[1]=='plassembler':
    assert args[-2:]==['--args','--no_chromosome'], args
    assert os.environ['PLASSEMBLER_DB'].endswith('/database')
else:
    assert '--args' not in args and '--no_chromosome' not in args, args
if args[1]=='flye': sys.exit(1)
Path(args[args.index('--out_prefix')+1]+'.fasta').write_text('>contig\nACGT\n')
""")
            executable.chmod(0o755)
            command=SOURCE.split('rule autocycler_candidate_assembly:')[1].split('"""')[1]
            command=command.replace('./scripts/',str(ROOT/'workflow/scripts')+'/')
            quote=lambda value: shlex.quote(str(value))
            env=dict(os.environ,PATH=str(p)+os.pathsep+os.environ['PATH'])
            jobs=[]
            for assembler in ['canu','flye','metamdbg','miniasm','necat','nextdenovo','plassembler','raven']:
                output=p/'jobs'/assembler/'sample_01'
                jobs.append(output)
                rendered=command.replace(':q}', '}').format(
                    input=SimpleNamespace(genome_size=quote(p/'size')),
                    output=quote(output),params=SimpleNamespace(read_type='ont_r10',reads=quote(subsets/'sample_01.fastq'),plassembler_db=quote(p/'database')),
                    wildcards=SimpleNamespace(candidate_assembler=assembler),threads=2,log=quote(p/f'{assembler}.log'))
                result=subprocess.run(['bash','-euo','pipefail','-c',rendered],env=env,capture_output=True)
                if assembler == 'flye':
                    self.assertEqual(result.returncode,0,result.stderr)
                    self.assertEqual((output/'status.txt').read_text().strip(),'failed')
                else:
                    self.assertEqual(result.returncode,0,result.stderr)
                    self.assertEqual((output/'status.txt').read_text().strip(),'success')
            import sys
            sys.path.insert(0,str(ROOT/'workflow/scripts'))
            from collect_autocycler_candidates import collect
            collect(jobs,p/'candidates')
            self.assertEqual(len(list((p/'candidates/assemblies').glob('*.fasta'))),7)
            failed=p/'jobs/failed/sample_01'
            failed.mkdir(parents=True)
            (failed/'status.txt').write_text('failed\n')
            collect([failed],p/'empty_candidates',allow_insufficient=True)
            self.assertEqual(list((p/'empty_candidates/assemblies').glob('*.fasta')),[])
            with self.assertRaisesRegex(ValueError,'multiple successful'):
                collect([failed],p/'too_few')
