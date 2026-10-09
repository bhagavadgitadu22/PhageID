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
            command=SOURCE.split('rule autocycler_candidate_assembly:')[1].split('"""')[1]
            quote=lambda value: shlex.quote(str(value))
            env=dict(os.environ,PATH=str(p)+os.pathsep+os.environ['PATH'])
            jobs=[]
            for assembler in ['canu','flye','metamdbg','miniasm','necat','nextdenovo','plassembler','raven']:
                output=p/'jobs'/assembler/'sample_01'
                jobs.append(output)
                rendered=command.replace(':q}', '}').format(
                    input=SimpleNamespace(genome_size=quote(p/'size')),
                    output=quote(output),params=SimpleNamespace(read_type='ont_r10',reads=quote(subsets/'sample_01.fastq')),
                    wildcards=SimpleNamespace(candidate_assembler=assembler),threads=2,log=quote(p/f'{assembler}.log'))
                subprocess.run(['bash','-euo','pipefail','-c',rendered],env=env,check=True,capture_output=True)
                self.assertEqual((output/'status.txt').read_text().strip(),'failed' if assembler=='flye' else 'success')
            import sys
            sys.path.insert(0,str(ROOT/'workflow/scripts'))
            from collect_autocycler_candidates import collect
            collect(jobs,p/'candidates')
            self.assertEqual(len(list((p/'candidates/assemblies').glob('*.fasta'))),7)
            self.assertIn('failed',(p/'candidates/candidate_status.tsv').read_text())
            with self.assertRaisesRegex(ValueError,'multiple successful'):
                collect([jobs[1]],p/'too_few')
