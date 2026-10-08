import os
from pathlib import Path
from types import SimpleNamespace
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class OptionalHostTests(unittest.TestCase):
    def test_sample_sheet_accepts_missing_host(self):
        source = (ROOT / 'workflow/Snakefile').read_text()
        source = source[source.index('with open(SAMPLES_FILE)'):source.index('def active_viral_samples')]
        with tempfile.TemporaryDirectory() as directory:
            sheet = Path(directory) / 'samples.tsv'
            sheet.write_text('a\t/reads/a\nb\t/reads/b\t\nc\t/reads/c\tNA\nd\t/reads/d\tNone\ne\t/reads/e\t-\nf\t/reads/f\t/host.fa\n')
            namespace = {'SAMPLES_FILE': str(sheet)}
            exec(source, namespace)
            self.assertEqual(namespace['HOSTS_LIST'], dict(a='', b='', c='', d='', e='', f='/host.fa'))
            self.assertEqual(len(namespace['READ_FILES']), 6)

    def test_missing_host_preserves_reads_and_logs_skip(self):
        source = (ROOT / 'workflow/rules/1_preprocessing_reads.smk').read_text()
        source = source[source.index('rule remove_bacterial_contamination:'):]
        command = source.split('"""')[1]
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            reads, output, log = (str(folder / name) for name in ['reads.fastq', 'cleaned.fastq', 'filter.log'])
            Path(reads).write_text('@read\nACGT\n+\nIIII\n')
            # Empty inputs render to an empty string in Snakemake shell formatting.
            command = command.replace(':q}', '}').format(
                params=SimpleNamespace(has_host=0), input=SimpleNamespace(ref='', reads=reads),
                output=output, log=log, threads=1)
            subprocess.run(['bash', '-euo', 'pipefail', '-c', command], check=True)
            self.assertEqual(Path(output).read_bytes(), Path(reads).read_bytes())
            self.assertIn('No host genome supplied', Path(log).read_text())


if __name__ == '__main__':
    unittest.main()
