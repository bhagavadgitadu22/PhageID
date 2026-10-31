import gzip
import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / 'workflow/Snakefile').read_text()
HELPER = SOURCE[SOURCE.index('def read_input_files'):SOURCE.index('def sample_reads')]
namespace = {'os': os}
exec(HELPER, namespace)
read_input_files = namespace['read_input_files']
namespace = {}
RULE_SOURCE = (ROOT / 'workflow/rules/01_preprocessing_reads.smk').read_text()
CONCAT = RULE_SOURCE.split('rule concat_reads:')[1].split('"""')[1]


class ReadInputTests(unittest.TestCase):
    def combine(self, path, output):
        files = read_input_files(str(path))
        command = CONCAT.replace(':q}', '}').format(
            input=shlex.join(files), output=shlex.quote(str(output)),
            log=shlex.quote(str(output) + '.log'))
        return subprocess.run(['bash', '-euo', 'pipefail', '-c', command], capture_output=True)

    def test_direct_plain_and_gzip_files(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            reads = b'@read\nACGT\n+\nIIII\n'
            for name in ['one read.fastq', 'one read.fq.gz']:
                source = p / name
                if name.endswith('.gz'):
                    with gzip.open(source, 'wb') as handle:
                        handle.write(reads)
                else:
                    source.write_bytes(reads)
                self.assertEqual(self.combine(source, p / 'output.fastq').returncode, 0)
                self.assertEqual((p / 'output.fastq').read_bytes(), reads)

    def test_directory_combines_only_fastq_files_in_order(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            directory = p / 'read directory'; directory.mkdir()
            first = b'@a\nACGT\n+\nIIII\n'
            second = b'@b\nTTTT\n+\nIIII\n'
            (directory / '02 reads.fq').write_bytes(second)
            with gzip.open(directory / '01 reads.fastq.gz', 'wb') as handle:
                handle.write(first)
            (directory / 'notes.txt').write_text('ignore this')
            (directory / 'nested.fastq').mkdir()
            self.assertEqual(self.combine(directory, p / 'output.fastq').returncode, 0)
            self.assertEqual((p / 'output.fastq').read_bytes(), first + second)

    def test_empty_directory_and_bad_gzip_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            with self.assertRaisesRegex(ValueError, 'No FASTQ files'):
                read_input_files(str(p))
            broken = p / 'broken.fastq.gz'; broken.write_bytes(b'not gzip')
            self.assertNotEqual(self.combine(broken, p / 'output.fastq').returncode, 0)
            empty = p / 'empty.fastq'; empty.touch()
            self.assertNotEqual(self.combine(empty, p / 'output.fastq').returncode, 0)

    def test_preflight_reports_all_missing_inputs(self):
        start = SOURCE.index('def validate_inputs():')
        end = SOURCE.index('\nvalidate_inputs()\n', start)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            empty_reads = root / 'empty_reads'
            empty_reads.mkdir()
            database_keys = ('genomad_db', 'checkv_db', 'pharokka_db', 'phold_db',
                             'plassembler_db', 'checkamg_db', 'vcontact3_db')
            config = {key: str(root / key) for key in database_keys}
            for key in database_keys:
                if key != 'checkv_db':
                    (root / key).mkdir()
            scope = dict(os=os, PHAGES_LIST=['sample'],
                         READ_FILES={'sample': str(root / 'missing.fastq')},
                         SHORT_READS={'sample': str(empty_reads)},
                         HOSTS_LIST={'sample': str(root / 'missing_host.fasta')},
                         config=config, read_input_files=read_input_files)
            exec(SOURCE[start:end], scope)
            with self.assertRaisesRegex(ValueError, 'Missing pipeline inputs') as error:
                scope['validate_inputs']()
            message = str(error.exception)
            for expected in ('sample long reads', 'sample short reads',
                             'sample host reference', 'checkv_db',
                             'genomad_marker_metadata.tsv'):
                self.assertIn(expected, message)


if __name__ == '__main__':
    unittest.main()
