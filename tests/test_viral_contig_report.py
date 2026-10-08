import csv
import importlib.util
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('report', ROOT/'workflow/scripts/viral_contig_report.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_csv(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class ViralReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name)
        (self.p/'flye.tsv').write_text('#seq_name\tlength\tcov.\tcirc.\trepeat\ncontig_1\t300\t42\tY\tN\ncontig_2\t80\t17\tN\tN\n')
        (self.p/'input.fa').write_text('>sample_with_underscores_contig_1\n'+'A'*90+'\n>sample_with_underscores_contig_2\n'+'C'*80+'\n')
        self.concat=[dict(contig_id='sample_with_underscores_contig_1',original_length=300,corrected_length=100,num_copies=3,status='corrected'),dict(contig_id='sample_with_underscores_contig_2',original_length=80,corrected_length=80,num_copies=0,status='no_repeats')]
        self.dtr=[dict(contig_id='sample_with_underscores_contig_1',original_length=100,corrected_length=90),dict(contig_id='sample_with_underscores_contig_2',original_length=80,corrected_length=80)]
        write_csv(self.p/'concat.csv',self.concat,list(self.concat[0]))
        write_csv(self.p/'dtr.csv',self.dtr,list(self.dtr[0]))

    def report(self):
        module.sample_report('sample_with_underscores',self.p/'input.fa',self.p/'flye.tsv',self.p/'concat.csv',self.p/'dtr.csv',self.p/'stats.tsv')

    def test_metadata_original_and_corrected_lengths(self):
        self.report()
        rows=read_tsv(self.p/'stats.tsv')
        self.assertEqual([row['Total bp'] for row in rows],['300','80'])
        self.assertEqual([row['Total corrected bp'] for row in rows],['90','80'])
        self.assertEqual(rows[0]['Coverage'],'42')
        self.assertEqual(rows[0]['Circular'],'Y')
        self.assertEqual(rows[0]['Number of concatemers broken'],'2')
        self.assertEqual(rows[0]['DTR length removed'],'10')
        self.assertEqual(rows[1]['DTR length removed'],'0')

    def test_summary_keeps_empty_samples_and_sums_both_lengths(self):
        self.report()
        module.write_tsv(self.p/'empty.tsv',module.CONTIG_FIELDS,[])
        module.combine_reports([('sample_with_underscores',self.p/'stats.tsv'),('empty',self.p/'empty.tsv')],self.p/'global')
        rows=read_tsv(self.p/'global/samples.tsv')
        self.assertEqual((rows[0]['Total bp'],rows[0]['Total corrected bp']),('380','170'))
        self.assertEqual(rows[1],dict(zip(module.SAMPLE_FIELDS,['empty','0','0','0','no_viral_contigs'])))
        self.assertEqual(len(read_tsv(self.p/'global/viral_contigs.tsv')),2)

    def test_empty_corrected_fasta_has_header_only_report(self):
        (self.p/'input.fa').write_text('')
        write_csv(self.p/'concat.csv',[],list(self.concat[0]))
        write_csv(self.p/'dtr.csv',[],list(self.dtr[0]))
        self.report()
        self.assertEqual(read_tsv(self.p/'stats.tsv'),[])

    def test_missing_metadata_fails(self):
        (self.p/'flye.tsv').write_text('#seq_name\tlength\tcov.\tcirc.\n')
        with self.assertRaisesRegex(ValueError,'No Flye metadata'): self.report()

    def test_checkpoint_filters_targets_and_skips_all_empty_comparison(self):
        module.write_tsv(self.p/'samples.tsv',module.SAMPLE_FIELDS,[dict(zip(module.SAMPLE_FIELDS,['ready','1','10','10','ready'])),dict(zip(module.SAMPLE_FIELDS,['empty','0','0','0','no_viral_contigs']))])
        checkpoint=SimpleNamespace(get=lambda: SimpleNamespace(output=SimpleNamespace(report=str(self.p))))
        def expand(patterns, sample):
            patterns=[patterns] if isinstance(patterns,str) else patterns
            samples=[sample] if isinstance(sample,str) else sample
            return [pattern.format(sample=name) for name in samples for pattern in patterns]
        namespace={'os':os,'checkpoints':SimpleNamespace(viral_sample_report=checkpoint),'PHAGES_LIST':['ready','empty'],'expand':expand}
        source=(ROOT/'workflow/Snakefile').read_text()
        exec(source[source.index('def active_viral_samples'):source.index('### Rules to include')],namespace)
        self.assertEqual(namespace['active_sample_files']('{sample}/annotation')(None),['ready/annotation'])
        self.assertEqual(namespace['per_sample_pipeline_inputs'](SimpleNamespace(sample='empty')),[str(self.p)])
        module.write_tsv(self.p/'samples.tsv',module.SAMPLE_FIELDS,[dict(zip(module.SAMPLE_FIELDS,['empty','0','0','0','no_viral_contigs']))])
        self.assertEqual(namespace['comparison_pipeline_inputs'](None),[str(self.p)])


if __name__ == '__main__': unittest.main()
