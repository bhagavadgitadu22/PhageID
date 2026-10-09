import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


class AssemblySelectionTests(unittest.TestCase):
    def test_filtered_first_and_all_reads_rescue(self):
        source = (ROOT / 'workflow/Snakefile').read_text()
        source = source[source.index('def candidate_assembly'):source.index('def sample_mapper')]
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)
            outputs = lambda **values: SimpleNamespace(output=SimpleNamespace(items=lambda: values.items()))
            pattern = '{sample}/{readset}/{assembler}'
            rules = SimpleNamespace(
                fix_circular_viral_contigs_per_sample_candidate=outputs(corrected=pattern+'/corrected.fa'),
                genomad_candidate=outputs(summary=pattern+'/summary.tsv'),
                checkv_candidate=outputs(checkv_quality=pattern+'/quality.tsv'),
                assembly_reads_flye=SimpleNamespace(output=SimpleNamespace(info=pattern+'/assembly_info.txt')),
                subsample_reads_flye=SimpleNamespace(output=SimpleNamespace(reads='{sample}/flye.{readset}.fastq')),
                remove_bacterial_contamination=SimpleNamespace(output=['{sample}/cleaned.fastq']))
            filter_report = p / 'filter.tsv'
            def status(**kwargs):
                return SimpleNamespace(output=[str(filter_report)])
            def quality(sample, readset, assembler):
                return SimpleNamespace(output=[str(p / f'{readset}_{assembler}.tsv')])
            namespace = dict(os=os, RESULTS_DIR=folder, READ_FILES={'long': '/reads', 'short': ''},
                HOSTS_LIST={'long': '/host', 'short': '/host'}, rules=rules,
                sample_reads=lambda wc: f'{wc.sample}/all.fastq',
                checkpoints=SimpleNamespace(assembly_read_status=SimpleNamespace(get=status),
                                            candidate_viral_quality=SimpleNamespace(get=quality)))
            exec(source, namespace)
            cases = [
                ('long', 'High-quality', None, None, None, 5, 'filtered', 'flye'),
                ('long', 'Complete', None, None, None, 5, 'filtered', 'flye'),
                ('long', 'Low-quality', 'High-quality', None, None, 5, 'filtered', 'autocycler'),
                ('long', 'Low-quality', None, None, None, 5, 'filtered', 'flye'),
                ('long', None, None, 'High-quality', None, 5, 'all', 'flye'),
                ('long', None, None, 'Low-quality', 'High-quality', 5, 'all', 'autocycler'),
                ('long', None, None, None, None, 5, 'all', 'flye'),
                ('long', None, None, 'High-quality', None, 0, 'all', 'flye'),
                ('short', 'Low-quality', None, None, None, 5, 'filtered', 'spades'),
                ('short', None, None, 'Low-quality', None, 5, 'all', 'spades')]
            for sample, primary, autocycler, full_primary, full_auto, filtered_count, readset, assembler in cases:
                filter_report.write_text(f'All reads number\tBacterial free reads\n10\t{filtered_count}\n')
                primary_assembler = 'flye' if sample == 'long' else 'spades'
                for name, category in [(f'filtered_{primary_assembler}', primary), ('filtered_autocycler', autocycler),
                                       (f'all_{primary_assembler}', full_primary), ('all_autocycler', full_auto)]:
                    (p / f'{name}.tsv').write_text('contig_id\tcheckv_quality\n'+(f'c\t{category}\n' if category else ''))
                with self.subTest(sample=sample, primary=primary, full=full_primary, filtered_count=filtered_count):
                    selected = namespace['selected_candidate_inputs'](SimpleNamespace(sample=sample))
                    self.assertEqual(selected['corrected'], f'{sample}/{readset}/{assembler}/corrected.fa')
                    expected_reads = (f'{sample}/flye.{readset}.fastq' if assembler == 'flye' else
                                      f'{sample}/cleaned.fastq' if readset == 'filtered' else f'{sample}/all.fastq')
                    self.assertEqual(selected['used_reads'], expected_reads)

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
