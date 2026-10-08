import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'workflow/scripts/dereplication_report.py'


class DereplicationReportTests(unittest.TestCase):
    def run_report(self, folder):
        return subprocess.run([sys.executable, '-S', str(SCRIPT), '--fasta', str(folder/'all.fa'),
            '--clusters', str(folder/'clusters.tsv'), '--representatives', str(folder/'selected.fa'),
            '--ani', str(folder/'ani.tsv'), '--output', str(folder/'report.tsv')], capture_output=True, text=True)

    def test_members_lengths_selection_and_singletons(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'all.fa').write_text('>short\nAC\n>long\nACGT\n>singleton\nAAA\n')
            (folder/'selected.fa').write_text('>long\nACGT\n>singleton\nAAA\n')
            (folder/'clusters.tsv').write_text('long\tlong,short\nsingleton\tsingleton\n')
            (folder/'ani.tsv').write_text('qname\ttname\tnum_alns\tpid\tqcov\ttcov\nlong\tshort\t1\t97.5\t50\t100\nshort\tlong\t1\t96\t90\t45\n')
            result = self.run_report(folder)
            self.assertEqual(result.returncode, 0, result.stderr)
            with open(folder/'report.tsv') as handle:
                rows = list(csv.reader(handle, delimiter='\t'))
            self.assertEqual(rows[0], ['Viral contig','Length','Selected','Representative','Representative length','Cluster size','ANI identity (%)','Aligned fraction (%)','Representative aligned fraction (%)'])
            self.assertEqual(rows[1:], [['short','2','FALSE','long','4','2','97.5','100','50'],['long','4','TRUE','long','4','2','100','100','100'],['singleton','3','TRUE','singleton','3','1','100','100','100']])
            (folder/'ani.tsv').write_text('qname\ttname\tnum_alns\tpid\tqcov\ttcov\nshort\tlong\t1\t96\t90\t45\n')
            self.assertIn('Missing representative-to-contig ANI', self.run_report(folder).stderr)
            (folder/'clusters.tsv').write_text('long\tlong,short,short\nsingleton\tsingleton\n')
            self.assertNotEqual(self.run_report(folder).returncode, 0)

    def test_empty_inputs_write_header(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for name in ['all.fa','selected.fa','clusters.tsv','ani.tsv']:
                (folder/name).write_text('')
            result = self.run_report(folder)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len((folder/'report.tsv').read_text().splitlines()), 1)
