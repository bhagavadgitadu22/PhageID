import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / 'workflow/scripts/plot_long_viral_reads.R'


class ViralReadPlotTests(unittest.TestCase):
    def test_empty_and_sparse_summaries_produce_distinct_pngs(self):
        rscript = os.environ.get('PHAGEID_RSCRIPT') or shutil.which('Rscript')
        if not rscript:
            self.skipTest('Rscript is not installed')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            images = []
            for name, reads in (('empty', ''),
                                ('sparse', 'read_1\t63000\nread_2\t63100\n')):
                summary = root / f'{name}.tsv'
                image = root / f'{name}.png'
                summary.write_text('seq_name\tlength\n' + reads)
                result = subprocess.run(
                    [rscript, str(SCRIPT), str(summary), str(image)],
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(image.read_bytes().startswith(b'\x89PNG\r\n\x1a\n'))
                images.append(image.read_bytes())
            self.assertNotEqual(images[0], images[1])


if __name__ == '__main__':
    unittest.main()
