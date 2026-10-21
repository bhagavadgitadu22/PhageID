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
        (self.p/'flye.tsv').write_text('#seq_name\tlength\tcov.\tcov_without_bacteria\tflye_circularity\nsample_with_underscores_contig_1\t300\t42\t30\tY\nsample_with_underscores_contig_2\t80\t17\t10\tN\n')
        (self.p/'input.fa').write_text('>sample_with_underscores_contig_1\n'+'A'*90+'\n>sample_with_underscores_contig_2\n'+'C'*80+'\n')
        self.concat=[dict(contig_id='sample_with_underscores_contig_1',original_length=300,corrected_length=100,num_copies=3,status='corrected'),dict(contig_id='sample_with_underscores_contig_2',original_length=80,corrected_length=80,num_copies=0,status='no_repeats')]
        self.dtr=[dict(contig_id='sample_with_underscores_contig_1',original_length=100,corrected_length=90),dict(contig_id='sample_with_underscores_contig_2',original_length=80,corrected_length=80)]
        write_csv(self.p/'concat.csv',self.concat,list(self.concat[0]))
        write_csv(self.p/'dtr.csv',self.dtr,list(self.dtr[0]))

    def concatenate_reports(self, reports, output_dir):
        import shlex
        import subprocess
        source='\n'.join(path.read_text() for path in sorted((ROOT/'workflow/rules').glob('*.smk')))
        command=source.split('checkpoint viral_report_global:')[1].split('"""')[1]
        quote=lambda value: shlex.quote(str(value))
        command=command.replace(':q}', '}').format(
            input=SimpleNamespace(contigs=shlex.join([str(report[1]) for report in reports]),
                                  summaries=shlex.join([str(report[2]) for report in reports])),
            output=SimpleNamespace(contigs=quote(output_dir/'assembly_report_per_viral_contig.tsv'),
                                   samples=quote(output_dir/'assembly_report_per_sample.tsv'),
                                   attempts=quote(output_dir/'assembly_attempts.tsv')),
            params=SimpleNamespace(attempts_script=quote(ROOT/'workflow/scripts/assembly_attempt_report.py'), results_dir=quote(self.p/'results'), contig_header=quote('\t'.join(module.CONTIG_FIELDS)),
                                   sample_header=quote('\t'.join(module.SAMPLE_FIELDS))),
            log=quote(self.p/'global.log'))
        subprocess.run(['bash','-euo','pipefail','-c',command],check=True)

    def report(self):
        module.sample_report('sample_with_underscores',self.p/'input.fa',self.p/'flye.tsv',self.p/'concat.csv',self.p/'dtr.csv',self.p/'stats.tsv', summary_output=self.p/'summary.tsv')

    def test_read_counts_in_both_reports_including_empty_samples(self):
        stats=self.p/'read_stats.tsv'
        filtering=self.p/'filter.tsv'
        stats.write_text('Assembly reads\tPercentage bacterial reads\tAll reads number\tReads used for assembly number\nbacterial free\t20\t100\t30\n')
        filtering.write_text('All reads number\tBacterial free reads\tPercentage bacterial reads\n100\t80\t20\n')
        for empty in [False, True]:
            if empty:
                (self.p/'input.fa').write_text('')
                write_csv(self.p/'concat.csv',[],list(self.concat[0]))
                write_csv(self.p/'dtr.csv',[],list(self.dtr[0]))
            module.sample_report('sample_with_underscores',self.p/'input.fa',self.p/'flye.tsv',self.p/'concat.csv',self.p/'dtr.csv',self.p/'stats.tsv',summary_output=self.p/'summary.tsv',assembler='flye',read_stats=stats,filter_report=filtering)
            summary=read_tsv(self.p/'summary.tsv')[0]
            self.assertEqual(list(summary)[1:6], ['Number_reads','Number_reads_without_host','Assembler','Assembly_method','Assembly_used_reads'])
            self.assertEqual([summary[key] for key in list(summary)[1:6]], ['100','80','flye','bacterial free','30'])
            for row in read_tsv(self.p/'stats.tsv'):
                self.assertEqual(list(row)[2:7], ['Number_reads','Number_reads_without_host','Assembler','Assembly_method','Assembly_used_reads'])
                self.assertEqual([row[key] for key in list(row)[2:7]], ['100','80','flye','bacterial free','30'])
                self.assertNotIn('Percentage bacterial reads', row)

    def test_failed_assembly_summary(self):
        row = module.sample_summary('failed', [], 'failed')
        self.assertEqual(row['Status'], 'assembly_failed')
        self.assertEqual(row['Viral_contigs'], 0)
        self.assertEqual(module.sample_summary('empty', [])['Status'], 'no_viral_contigs')

    def test_metadata_original_and_corrected_lengths(self):
        self.report()
        rows=read_tsv(self.p/'stats.tsv')
        self.assertEqual([row['Total_bp'] for row in rows],['300','80'])
        self.assertEqual([row['Total_corrected_bp'] for row in rows],['90','80'])
        self.assertEqual(rows[0]['Assembler'],'flye')
        self.assertEqual(rows[0]['Coverage'],'42')
        self.assertEqual(rows[0]['Coverage_without_bacteria'],'30')
        self.assertEqual(rows[0]['Flye_circularity'],'Y')
        self.assertEqual(rows[0]['Number_of_concatemers_broken'],'2')
        self.assertEqual(rows[0]['DTR_length_removed'],'10')
        self.assertEqual(rows[1]['DTR_length_removed'],'0')

    def test_summary_keeps_empty_samples_and_sums_both_lengths(self):
        self.report()
        module.write_tsv(self.p/'empty.tsv',module.CONTIG_FIELDS,[])
        module.write_tsv(self.p/'empty_summary.tsv',module.SAMPLE_FIELDS,[module.sample_summary('empty',[])])
        self.concatenate_reports([('sample_with_underscores',self.p/'stats.tsv',self.p/'summary.tsv'),('empty',self.p/'empty.tsv',self.p/'empty_summary.tsv')],self.p/'global')
        rows=read_tsv(self.p/'global/assembly_report_per_sample.tsv')
        self.assertEqual((rows[0]['Total_bp'],rows[0]['Total_corrected_bp']),('380','170'))
        self.assertEqual(rows[1],dict(zip(module.SAMPLE_FIELDS,['empty','NA','NA','NA','NA','NA','0','0','0','no_viral_contigs'])))
        self.assertEqual(len(read_tsv(self.p/'global/assembly_report_per_viral_contig.tsv')),2)

    def test_manual_length_change_keeps_dtr_count_separate(self):
        (self.p/'input.fa').write_text('>sample_with_underscores_contig_1\n'+'A'*75+'\n>sample_with_underscores_contig_2\n'+'C'*80+'\n')
        manual=self.p/'manual.tsv'
        module.write_tsv(manual, ['Viral_contig', 'Auto_corrected_bp', 'Final_bp',
                                  'Manual_correction', 'Manual_length_change_bp'], [
            dict(Viral_contig='sample_with_underscores_contig_1', Auto_corrected_bp=90,
                 Final_bp=75, Manual_correction='Yes', Manual_length_change_bp=-15),
            dict(Viral_contig='sample_with_underscores_contig_2', Auto_corrected_bp=80,
                 Final_bp=80, Manual_correction='No', Manual_length_change_bp=0)])
        module.sample_report('sample_with_underscores',self.p/'input.fa',self.p/'flye.tsv',
                             self.p/'concat.csv',self.p/'dtr.csv',self.p/'stats.tsv',
                             manual_report=manual)
        row=read_tsv(self.p/'stats.tsv')[0]
        self.assertEqual((row['Total_corrected_bp'],row['DTR_length_removed'],
                          row['Manual_correction'],row['Manual_length_change_bp']),
                         ('75','10','Yes','-15'))

    def test_empty_corrected_fasta_has_header_only_report(self):
        (self.p/'input.fa').write_text('')
        write_csv(self.p/'concat.csv',[],list(self.concat[0]))
        write_csv(self.p/'dtr.csv',[],list(self.dtr[0]))
        self.report()
        self.assertEqual(read_tsv(self.p/'stats.tsv'),[])

    def test_enrichment_joins_checkv_and_genomad_ids(self):
        self.report()
        columns = ['contig_id', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'miuvig_quality', 'completeness', 'contamination']
        ids = [row['contig_id'] for row in self.concat]
        module.write_tsv(self.p/'checkv.tsv', columns, [dict(zip(columns,[name,'10','8','2','High-quality','High-quality','95.2','3.4'])) for name in ids])
        module.write_tsv(self.p/'genomad.tsv', ['seq_name','topology','taxonomy'], [dict(seq_name='contig_1',topology='Provirus',taxonomy='Viruses;Caudoviricetes'),dict(seq_name='contig_2',topology='DTR',taxonomy='Unclassified')])
        module.sample_report('sample_with_underscores',self.p/'input.fa',self.p/'flye.tsv',self.p/'concat.csv',self.p/'dtr.csv',self.p/'enriched.tsv',self.p/'checkv.tsv',self.p/'genomad.tsv',self.p/'enriched_summary.tsv')
        rows=read_tsv(self.p/'enriched.tsv')
        self.assertEqual([row['geNomad_provirus'] for row in rows], ['Yes','No'])
        self.assertEqual(rows[0]['CheckV_completeness'], '95.2')
        self.assertEqual(rows[0]['CheckV_contamination'], '3.4')
        fields = list(rows[0])
        self.assertEqual(fields[fields.index('CheckV_completeness') + 1], 'CheckV_contamination')
        self.assertNotIn('CheckV completeness_method', fields)
        self.assertEqual(rows[0]['CheckV_gene_count'], '10')
        self.assertEqual(rows[0]['geNomad_taxonomy'], 'Viruses;Caudoviricetes')
        self.concatenate_reports([('sample_with_underscores',self.p/'enriched.tsv',self.p/'enriched_summary.tsv')],self.p/'enriched_global')
        self.assertEqual(read_tsv(self.p/'enriched_global/assembly_report_per_viral_contig.tsv'), rows)

    def test_provirus_uses_parent_coverage_and_region_length(self):
        name='sample_with_underscores_contig_1|provirus_11_110'
        (self.p/'input.fa').write_text('>'+name+'\n'+'A'*100+'\n')
        concat=dict(contig_id=name,original_length=100,corrected_length=100,num_copies=0,status='no_repeats')
        dtr=dict(contig_id=name,original_length=100,corrected_length=100)
        write_csv(self.p/'concat.csv',[concat],list(concat))
        write_csv(self.p/'dtr.csv',[dtr],list(dtr))
        (self.p/'flye.tsv').write_text('#seq_name\tlength\tcov.\tcov_without_bacteria\tflye_circularity\n'+name+'\t100\t42\t30\tNA\n')
        self.report()
        row=read_tsv(self.p/'stats.tsv')[0]
        self.assertEqual((row['Total_bp'],row['Coverage'],row['Flye_circularity']),('100','42','NA'))

    def test_empty_enrichment_needs_no_annotations(self):
        module.write_tsv(self.p/'empty.tsv',module.CONTIG_FIELDS,[])
        (self.p/'input.fa').write_text('')
        write_csv(self.p/'concat.csv',[],list(self.concat[0]))
        write_csv(self.p/'dtr.csv',[],list(self.dtr[0]))
        module.sample_report('empty',self.p/'input.fa',self.p/'flye.tsv',self.p/'concat.csv',self.p/'dtr.csv',self.p/'enriched.tsv')
        self.assertEqual(read_tsv(self.p/'enriched.tsv'), [])

    def test_per_sample_then_global_empty_reports(self):
        (self.p/'input.fa').write_text('')
        write_csv(self.p/'concat.csv',[],list(self.concat[0]))
        write_csv(self.p/'dtr.csv',[],list(self.dtr[0]))
        contigs=self.p/'results/empty/assembly_stats.tsv'
        summary=self.p/'results/empty/viral_sample_report.tsv'
        module.sample_report('empty',self.p/'input.fa',self.p/'flye.tsv',self.p/'concat.csv',self.p/'dtr.csv',contigs,summary_output=summary)
        self.concatenate_reports([('empty',contigs,summary)],self.p/'global')
        self.assertEqual(read_tsv(self.p/'results/empty/assembly_stats.tsv'), [])
        self.assertEqual(read_tsv(self.p/'global/assembly_report_per_viral_contig.tsv'), [])
        self.assertEqual(read_tsv(self.p/'global/assembly_report_per_sample.tsv')[0]['Status'], 'no_viral_contigs')

    def test_global_concatenation_with_no_samples_writes_headers(self):
        self.concatenate_reports([],self.p/'global')
        self.assertEqual(read_tsv(self.p/'global/assembly_report_per_viral_contig.tsv'), [])
        self.assertEqual(read_tsv(self.p/'global/assembly_report_per_sample.tsv'), [])

    def test_fasta_parser_multiline_and_duplicate_ids(self):
        (self.p/'lengths.fa').write_text('>one description\nACG\nTT\n>two\nC\n')
        self.assertEqual(module.fasta_lengths(self.p/'lengths.fa'),dict(one=5,two=1))
        (self.p/'lengths.fa').write_text('>one\nAC\n>one\nGT\n')
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            module.fasta_lengths(self.p/'lengths.fa')

    def test_checkv_skips_empty_fasta(self):
        import subprocess
        source='\n'.join(path.read_text() for path in sorted((ROOT/'workflow/rules').glob('*.smk')))
        command=source.split('rule checkv_candidate:')[1].split('"""')[1]
        (self.p/'empty.fa').write_text('')
        output=self.p/'checkv/quality_summary.tsv'
        stale_genes=self.p/'checkv/tmp/proteins.faa'
        stale_genes.parent.mkdir(parents=True)
        stale_genes.write_text('>old_contig_gene\nMKK\n')
        output.write_text('outdated result\n')
        log=self.p/'checkv.log'
        command=command.replace(':q}', '}').format(threads=1,input=SimpleNamespace(assembly=self.p/'empty.fa',db=self.p/'unused_db'),output=SimpleNamespace(checkv_quality=output),log=log)
        subprocess.run(['bash','-euo','pipefail','-c',command],check=True)
        self.assertEqual(read_tsv(output), [])
        self.assertFalse(stale_genes.exists())
        self.assertIn('skipping CheckV',log.read_text())

    def test_missing_metadata_fails(self):
        (self.p/'flye.tsv').write_text('#seq_name\tlength\tcov.\tcirc.\n')
        with self.assertRaisesRegex(ValueError,'No assembly metadata'): self.report()

    def test_checkpoint_filters_targets_and_skips_all_empty_comparison(self):
        sample_report=self.p/'assembly_report_per_sample.tsv'
        module.write_tsv(sample_report,module.SAMPLE_FIELDS,[dict(zip(module.SAMPLE_FIELDS,['ready','100','90','flye','bacterial free','80','1','10','10','ready'])),dict(zip(module.SAMPLE_FIELDS,['empty','NA','NA','NA','NA','NA','0','0','0','no_viral_contigs']))])
        checkpoint=SimpleNamespace(get=lambda: SimpleNamespace(output=SimpleNamespace(samples=str(sample_report))))
        def expand(patterns, sample):
            patterns=[patterns] if isinstance(patterns,str) else patterns
            samples=[sample] if isinstance(sample,str) else sample
            return [pattern.format(sample=name) for name in samples for pattern in patterns]
        namespace={'os':os,'checkpoints':SimpleNamespace(viral_report_global=checkpoint),'PHAGES_LIST':['ready','empty'],'expand':expand}
        source=(ROOT/'workflow/Snakefile').read_text()
        exec(source[source.index('def active_viral_samples'):source.index('### Rules to include')],namespace)
        self.assertEqual(namespace['active_sample_files']('{sample}/annotation')(None),['ready/annotation'])
        self.assertEqual(namespace['per_sample_pipeline_inputs'](SimpleNamespace(sample='empty')),[str(sample_report)])
        module.write_tsv(sample_report,module.SAMPLE_FIELDS,[dict(zip(module.SAMPLE_FIELDS,['empty','NA','NA','NA','NA','NA','0','0','0','no_viral_contigs']))])
        self.assertEqual(namespace['comparison_pipeline_inputs'](None),[str(sample_report)])


if __name__ == '__main__': unittest.main()
