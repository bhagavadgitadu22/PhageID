import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('subsample_flye', ROOT / 'workflow/scripts/subsample_flye_reads.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FlyeSubsamplingTests(unittest.TestCase):
    def test_sampling_is_reproducible_and_keeps_whole_records(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            reads = ''.join(f'@r{i}\n'+ 'A' * (i + 5) + '\n+\n' + 'I' * (i + 5) + '\n' for i in range(30))
            (p / 'input.fastq').write_text(reads)
            for name in ['first', 'second']:
                module.subsample(p / 'input.fastq', p / name, p / f'{name}.tsv', 10, 10, 42)
            self.assertEqual((p / 'first').read_bytes(), (p / 'second').read_bytes())
            records = list(module.fastq_records(p / 'first'))
            bases = sum(length for _, length in records)
            self.assertGreaterEqual(bases, 100)
            self.assertLess(bases, 100 + 34)
            self.assertLess(len(records), 30)
            self.assertEqual((p / 'input.fastq').read_text(), reads)
            with open(p / 'first.tsv') as handle:
                row = next(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(int(row['Selected bp']), bases)
            self.assertEqual(float(row['Estimated coverage']), bases / 10)

    def test_below_target_preserves_all_reads(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            reads = '@read\nACGT\n+\nIIII\n'
            (p / 'input').write_text(reads)
            module.subsample(p / 'input', p / 'output', p / 'report')
            self.assertEqual((p / 'output').read_text(), reads)

    def test_invalid_empty_inputs_and_settings_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            for contents in ['', '@read\nACGT\n+\nII\n']:
                (p / 'input').write_text(contents)
                with self.assertRaises(ValueError):
                    module.subsample(p / 'input', p / 'output', p / 'report')
            with self.assertRaisesRegex(ValueError, 'must be positive'):
                module.subsample(p / 'input', p / 'output', p / 'report', coverage=0)


if __name__ == '__main__':
    unittest.main()
