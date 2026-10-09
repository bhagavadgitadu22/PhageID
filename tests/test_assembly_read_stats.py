import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
from assembly_read_stats import filtering_stats, selected_stats


def row(path):
    with open(path) as handle:
        return next(csv.DictReader(handle, delimiter='\t'))


class ReadStatsTests(unittest.TestCase):
    def test_filtering_and_actual_subsample_counts(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            read = '@r\nACGT\n+\nIIII\n'
            (p / 'all').write_text(read * 10)
            (p / 'filtered').write_text(read * 6)
            (p / 'used').write_text(read * 2)
            filtering_stats(p / 'all', p / 'filtered', True, p / 'filter.tsv')
            self.assertEqual(float(row(p / 'filter.tsv')['Percentage bacterial reads']), 40)
            for readset, expected in [('filtered', 'bacterial free'), ('all', 'all')]:
                selected_stats(p / 'filter.tsv', p / 'used', readset, p / 'selected.tsv')
                report = row(p / 'selected.tsv')
                self.assertEqual(report['Assembly reads'], expected)
                self.assertEqual(report['All reads number'], '10')
                self.assertEqual(report['Reads used for assembly number'], '2')

    def test_no_host_has_unknown_bacterial_fraction_and_empty_filter_is_valid(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            (p / 'all').write_text('@r\nACGT\n+\nIIII\n')
            (p / 'filtered').write_text('')
            filtering_stats(p / 'all', p / 'all', False, p / 'filter.tsv')
            self.assertEqual(row(p / 'filter.tsv')['Percentage bacterial reads'], 'NA')
            filtering_stats(p / 'all', p / 'filtered', True, p / 'filter.tsv')
            self.assertEqual(row(p / 'filter.tsv')['Bacterial free reads'], '0')
            self.assertEqual(float(row(p / 'filter.tsv')['Percentage bacterial reads']), 100)
            with self.assertRaisesRegex(ValueError, 'assembly input reads'):
                selected_stats(p / 'filter.tsv', p / 'all', 'filtered', p / 'selected.tsv')


if __name__ == '__main__':
    unittest.main()
