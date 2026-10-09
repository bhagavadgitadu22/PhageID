import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


class AssemblySelectionTests(unittest.TestCase):
    def test_checkv_selects_route(self):
        source = (ROOT / 'workflow/Snakefile').read_text()
        source = source[source.index('def candidate_assembly'):source.index('def sample_mapper')]
        with tempfile.TemporaryDirectory() as folder:
            quality = Path(folder) / 'quality.tsv'
            outputs = lambda **values: SimpleNamespace(output=SimpleNamespace(items=lambda: values.items()))
            rules = SimpleNamespace(
                fix_circular_viral_contigs_per_sample_candidate=outputs(corrected='{sample}/{assembler}/corrected.fa'),
                genomad_candidate=outputs(summary='{sample}/{assembler}/summary.tsv'),
                checkv_candidate=outputs(checkv_quality='{sample}/{assembler}/quality.tsv'),
                assembly_reads_flye=SimpleNamespace(output=SimpleNamespace(info='{sample}/flye/assembly_info.txt')),
                mapped_assembly_info=SimpleNamespace(output=['{sample}/{assembler}/assembly_info.tsv']))
            checkpoint = SimpleNamespace(get=lambda **kwargs: SimpleNamespace(output=[str(quality)]))
            namespace = dict(os=os, RESULTS_DIR=folder, READ_FILES={'long': '/reads', 'short': ''},
                             rules=rules, checkpoints=SimpleNamespace(primary_viral_quality=checkpoint))
            exec(source, namespace)
            for category, expected in [('Complete', 'flye'), ('High-quality', 'flye'),
                                       ('Medium-quality', 'autocycler'), ('Low-quality', 'autocycler'),
                                       ('Not-determined', 'autocycler'), (None, 'autocycler')]:
                quality.write_text('contig_id\tcheckv_quality\n' +
                                   (f'a\t{category}\n' if category else ''))
                with self.subTest(category=category):
                    selected = namespace['selected_candidate_inputs'](SimpleNamespace(sample='long'))
                    self.assertEqual(selected['corrected'], f'long/{expected}/corrected.fa')
                    selected = namespace['selected_candidate_inputs'](SimpleNamespace(sample='short'))
                    self.assertEqual(selected['corrected'], 'short/spades/corrected.fa')
            quality.write_text('contig_id\tcheckv_quality\na\tLow-quality\nb\tHigh-quality\n')
            self.assertEqual(namespace['selected_candidate_inputs'](SimpleNamespace(sample='long'))['corrected'],
                             'long/flye/corrected.fa')

    def test_circular_coverage_counts_wrapped_bases_once(self):
        import pysam
        sys.path.insert(0, str(ROOT / 'workflow/scripts'))
        from assembly_info import bam_depths
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            for circular in [False, True]:
                bam_path = p / f'{circular}.bam'
                header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'a', 'LN': 100}, {'SN': 'zero', 'LN': 100}]}
                if circular:
                    header['CO'] = ['theBIGbam:circular=true']
                with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as bam:
                    for flag in [0, 256, 1024]:
                        read = pysam.AlignedSegment()
                        read.query_name = str(flag); read.query_sequence = 'A' * 20
                        read.flag = flag; read.reference_id = 0; read.reference_start = 90
                        read.cigar = [(0, 20)]; read.mapping_quality = 60
                        bam.write(read)
                depths = bam_depths(bam_path, {'a': 100, 'zero': 100})
                self.assertEqual(depths, {'a': 0.2 if circular else 0.1, 'zero': 0.0})
                with self.assertRaisesRegex(ValueError, 'BAM references'):
                    bam_depths(bam_path, {'a': 200})

    def test_short_read_sample_sheet(self):
        source = (ROOT / 'workflow/Snakefile').read_text()
        source = source[source.index('with open(SAMPLES_FILE)'):source.index('def active_viral_samples')]
        with tempfile.TemporaryDirectory() as folder:
            sheet = Path(folder) / 'samples.tsv'
            sheet.write_text('short\t\t/reads/single.fastq.gz\t/host.fa\nboth\t/long\t/short.fastq\tNA\n')
            namespace = {'SAMPLES_FILE': str(sheet)}
            exec(source, namespace)
            self.assertEqual(namespace['READ_FILES'], {'short': '', 'both': '/long'})
            self.assertEqual(namespace['SHORT_READS'], {'short': '/reads/single.fastq.gz', 'both': '/short.fastq'})
            self.assertEqual(namespace['HOSTS_LIST'], {'short': '/host.fa', 'both': ''})
            self.assertEqual(namespace['sample_mapper'](SimpleNamespace(sample='short')), 'minimap2-sr-secondary')
            self.assertEqual(namespace['sample_mapper'](SimpleNamespace(sample='both')), 'minimap2-ont')
            sheet.write_text('missing\t\t\t\n')
            with self.assertRaisesRegex(ValueError, 'No long or short reads'):
                exec(source, {'SAMPLES_FILE': str(sheet)})


if __name__ == '__main__':
    unittest.main()
