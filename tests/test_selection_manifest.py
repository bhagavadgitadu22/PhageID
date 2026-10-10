import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
import assembly_selection_manifest as selection


class SelectionManifestTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.sample = 'sample'
        self.candidate = self.root / 'sample/assembly_attempts/filtered/flye'
        self.sources = {key: self.candidate / value for key, value in selection.SOURCE_FILES.items()}
        self.sources.update(flye_info=self.candidate/'assembly_info.txt',
                            filter_report=self.root/'sample/reads/host_filtering_stats.tsv',
                            used_reads=self.root/'sample/reads/flye.sample.filtered.fastq')
        for path in self.sources.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        self.sources['assembly'].write_text('>c\nACGTACGT\n')
        self.sources['fasta'].write_text('>c\nACGTACGT\n')
        self.sources['corrected'].write_text('>sample_c\nACGTACGT\n')
        self.sources['status'].write_text('success\n')
        self.sources['used_reads'].write_text('@r\nACGT\n+\nIIII\n')
        self.sources['flye_info'].write_text('#seq_name\tlength\tcov.\tcirc.\nc\t8\t10\tY\n')
        self.sources['filter_report'].write_text('All reads number\tBacterial free reads\tPercentage bacterial reads\n2\t1\t50\n')
        for key in ['concatemer_report', 'dtr_report']:
            self.sources[key].write_text('contig_id,original_length,corrected_length\nsample_c,8,8\n')
        self.sources['checkv_quality'].write_text('contig_id\tcheckv_quality\nsample_c\tHigh-quality\n')
        self.sources['summary'].write_text('seq_name\ttopology\ttaxonomy\nc\tLinear\tViruses\n')
        self.manifest = self.root/'manifest.json'
        self.destinations = {key: self.root/'published'/key for key in selection.DESTINATION_KEYS}

    def test_create_and_publish_reuses_a_resolved_choice(self):
        selection.create(self.sample, list(self.sources.values()), self.manifest)
        manifest = json.loads(self.manifest.read_text())
        self.assertEqual((manifest['assembler'], manifest['readset']), ('flye', 'filtered'))
        for _ in range(2):
            selection.publish(self.sample, self.manifest, self.destinations)
            self.assertTrue(self.destinations['corrected'].read_text().startswith('>sample_c'))
            self.assertEqual(self.destinations['checkv_quality'].read_text(), self.sources['checkv_quality'].read_text())
            self.assertEqual(self.destinations['assembler'].read_text(), 'flye\n')
        with open(self.destinations['read_stats']) as handle:
            row = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(row['Reads used for assembly number'], '1')

    def test_unresolved_checkpoint_inputs_cannot_write_manifest(self):
        with self.assertRaisesRegex(ValueError, 'unresolved'):
            selection.create(self.sample, [self.sources['checkv_quality']]*11, self.manifest)
        self.assertFalse(self.manifest.exists())

    def test_tsv_in_place_of_fasta_cannot_overwrite_published_files(self):
        selection.create(self.sample, list(self.sources.values()), self.manifest)
        self.destinations['corrected'].parent.mkdir()
        self.destinations['corrected'].write_text('existing result')
        self.sources['corrected'].write_text(self.sources['checkv_quality'].read_text())
        with self.assertRaisesRegex(ValueError, 'FASTA header'):
            selection.publish(self.sample, self.manifest, self.destinations)
        self.assertEqual(self.destinations['corrected'].read_text(), 'existing result')
        self.assertFalse(self.destinations['assembler'].exists())

    def test_failed_attempt_remains_publishable(self):
        for key in ['assembly', 'fasta', 'corrected']:
            self.sources[key].write_text('')
        self.sources['status'].write_text('failed\n')
        self.sources['flye_info'].write_text('')
        for key in ['concatemer_report', 'dtr_report', 'checkv_quality', 'summary']:
            self.sources[key].write_text(self.sources[key].read_text().splitlines()[0]+'\n')
        selection.create(self.sample, list(self.sources.values()), self.manifest)
        selection.publish(self.sample, self.manifest, self.destinations)
        self.assertEqual(self.destinations['status'].read_text(), 'failed\n')
        self.assertEqual(self.destinations['corrected'].read_text(), '')

    def test_corrupt_reads_rejected_before_publication(self):
        selection.create(self.sample, list(self.sources.values()), self.manifest)
        self.sources['used_reads'].write_text('@r\nACGT\n+\nIIII\n@broken\nACGT\n')
        with self.assertRaises(ValueError):
            selection.publish(self.sample, self.manifest, self.destinations)
        self.assertFalse(self.destinations['corrected'].exists())


if __name__ == '__main__':
    unittest.main()
