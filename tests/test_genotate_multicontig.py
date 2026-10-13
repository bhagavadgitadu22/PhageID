"""Regression checks for multi-contig conversion and plot routing."""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'workflow' / 'scripts' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


converter = load('reformat_pharokka_genotate')
plotter = load('plotting_pharokka_genotate')
pharokka_plotter = load('plot_pharokka')


class MultiContigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / 'input.fa'
        self.fasta.write_text('>a\n' + 'ACGT' * 25 + '\n>b\n' + 'ACGT' * 50 + '\n>no_genes\nACGT\n')
        self.genes = self.root / 'genes.tsv'
        with self.genes.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['ID', 'length', 'phrog', 'annotation', 'category'])
            writer.writerow(['a_CDS_[1..9]', 9, 1, 'tail fiber; protein', 'tail'])
            writer.writerow(['b_CDS_[complement(10..18)]', 9, 2, 'capsid protein', 'head and packaging'])
        self.gff = self.root / 'out.gff'
        converter.convert(self.genes, self.fasta, self.gff)

    def test_converter_preserves_contigs_products_and_strands(self):
        features = plotter.read_features(self.gff, {'a': 100, 'b': 200, 'no_genes': 4})
        self.assertEqual(set(features), {'a', 'b'})
        self.assertEqual(features['a'][0].qualifiers['product'], ['tail fiber; protein'])
        self.assertEqual(features['a'][0].qualifiers['function'], ['tail'])
        self.assertEqual(features['b'][0].location.strand, -1)
        self.assertIn('##sequence-region no_genes 1 4', self.gff.read_text())

    def test_partial_boundaries_on_both_strands(self):
        cases = [('a_CDS_[<1..9]', '+', '.'),
                 ('a_CDS_[1..>9]', '+', '0'),
                 ('b_CDS_[complement(<1..18)]', '-', '0'),
                 ('b_CDS_[complement(10..>18)]', '-', '.'),
                 ('b_CDS_[complement(<1..>18)]', '-', '.')]
        with self.genes.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            for identifier, _, _ in cases:
                writer.writerow([identifier, 9, 1, 'capsid protein', 'head and packaging'])
        converter.convert(self.genes, self.fasta, self.gff)
        lines = [line.split('\t') for line in self.gff.read_text().split('##FASTA')[0].splitlines() if not line.startswith('#')]
        self.assertEqual(len(lines), len(cases))
        for fields, (identifier, strand, phase) in zip(lines, cases):
            self.assertEqual(fields[6:8], [strand, phase])
            self.assertIn('partial=true', fields[8])
            from urllib.parse import unquote
            self.assertIn('ID=' + identifier, unquote(fields[8]))
        features = plotter.read_features(self.gff, {'a': 100, 'b': 200, 'no_genes': 4})
        self.assertEqual(features['b'][0].location.start, 0)
        self.assertEqual(features['b'][0].location.end, 18)

    def test_partial_coordinates_still_require_valid_bounds_and_syntax(self):
        for identifier in ['a_CDS_[<0..9]', 'a_CDS_[1..>101]', 'a_CDS_[complement(<1..9]', 'a_CDS_[<1..9)]']:
            self.genes.write_text(identifier + '\t9\t1\tx\ttail\n')
            with self.subTest(identifier=identifier), self.assertRaises(ValueError):
                converter.convert(self.genes, self.fasta, self.gff)

    def test_pharokka_trailing_empty_columns_are_accepted(self):
        baseline = plotter.read_features(self.gff, {'a': 100, 'b': 200, 'no_genes': 4})
        lines = self.gff.read_text().splitlines()
        self.gff.write_text('\n'.join(line + '\t\t' if '\tCDS\t' in line else line for line in lines) + '\n')
        features = plotter.read_features(self.gff, {'a': 100, 'b': 200, 'no_genes': 4})
        self.assertEqual(set(features), set(baseline))
        for contig in features:
            self.assertEqual(features[contig][0].location, baseline[contig][0].location)
            self.assertEqual(features[contig][0].qualifiers, baseline[contig][0].qualifiers)

    def test_pharokka_trna_extra_anticodon_column_is_ignored(self):
        trna = ('a\tprofile:tRNAscan-SE:2.0.13\ttRNA\t10\t20\t39.0\t+\t.\t'
                'ID=tRNA_1;product=tRNA-Lys(AAA);anticodon=(pos:12..14,aa:Lys,seq:TTT)\tAAA\n')
        self.gff.write_text(self.gff.read_text().replace('##FASTA\n', trna + '##FASTA\n'))
        features = plotter.read_features(self.gff, {'a': 100, 'b': 200, 'no_genes': 4})
        self.assertEqual(len(features['a']), 1)

    def test_gff_extra_nonempty_columns_and_missing_columns_are_rejected(self):
        row = 'a\ttest\tCDS\t1\t9\t.\t+\t0\tID=one'
        for malformed in [row + '\textra', row.rsplit('\t', 1)[0]]:
            self.gff.write_text(malformed + '\n')
            with self.subTest(row=malformed), self.assertRaisesRegex(ValueError, 'nine GFF columns'):
                plotter.read_features(self.gff, {'a': 100})

    def test_unknown_contig_is_rejected(self):
        self.genes.write_text('missing_CDS_[1..3]\t3\t1\tx\ttail\n')
        with self.assertRaisesRegex(ValueError, 'unknown contig'):
            converter.convert(self.genes, self.fasta, self.gff)

    def test_plot_routes_each_contig_using_its_length(self):
        # Exercise routing without requiring the rendering dependency locally.
        calls, arrows = [], []
        class Track:
            def axis(self, **kwargs): pass
            def xticks_by_interval(self, **kwargs): pass
            def genomic_features(self, features, **kwargs): arrows.extend(features)
        class Sector:
            def text(self, *args, **kwargs): pass
            def add_track(self, *args): return Track()
        class Circos:
            def __init__(self, sectors): calls.append(sectors)
            def get_sector(self, name): return Sector()
            def plotfig(self):
                return SimpleNamespace(savefig=lambda path, **kwargs: Path(path).write_bytes(b'plot'))
        output = self.root / 'plots'
        args = ['plotter', '--fasta', str(self.fasta), '--genotate-gff', str(self.gff), '--pharokka-gff', str(self.gff), '--output-dir', str(output)]
        with patch.dict(sys.modules, {'pycirclize': SimpleNamespace(Circos=Circos)}), patch.object(sys, 'argv', args), patch.object(plotter.plt, 'close'):
            plotter.main()
        self.assertEqual(calls, [{'a': 100}, {'b': 200}])
        self.assertEqual(len(list(output.glob('*.png'))), 2)
        self.assertEqual(len(arrows), 4)  # Each CDS once in each annotation source.

    def test_genotate_plot_skips_unknown_only_contigs(self):
        self.genes.write_text('ID\tlength\tphrog\tannotation\tcategory\n'
                              'a_CDS_[1..9]\t9\tNo_PHROG\thypothetical protein\tunknown function\n')
        converter.convert(self.genes, self.fasta, self.gff)
        output = self.root / 'unknown_plots'
        output.mkdir()
        (output / 'stale_annotated_by_genotate_pharokka.png').write_bytes(b'old')
        args = ['plotter', '--fasta', str(self.fasta), '--genotate-gff', str(self.gff),
                '--pharokka-gff', str(self.gff), '--output-dir', str(output)]
        with patch.object(sys, 'argv', args):
            plotter.main()
        self.assertEqual(list(output.glob('*.png')), [])

    def test_pharokka_multiplot_filters_unlabelable_contigs(self):
        records = []
        for name, product in [('unknown', 'hypothetical protein'), ('known', 'tail protein')]:
            record = SeqRecord(Seq('ATG' * 10), id=name, name=name, description='')
            record.annotations['molecule_type'] = 'DNA'
            record.features = [SeqFeature(FeatureLocation(0, 9), type='CDS',
                                          qualifiers={'product': [product]})]
            records.append(record)
        genbank = self.root / 'pharokka.gbk'
        SeqIO.write(records, genbank, 'genbank')
        output = self.root / 'pharokka_plots'
        def fake_run(command, check):
            self.assertTrue(check)
            self.assertEqual([record.id for record in SeqIO.parse(command[3], 'genbank')], ['known'])
            Path(command[5]).mkdir()
            (Path(command[5]) / 'known.png').write_bytes(b'plot')
        with patch.object(pharokka_plotter.subprocess, 'run', side_effect=fake_run) as run:
            pharokka_plotter.plot(genbank, output)
            self.assertEqual(run.call_count, 1)
            self.assertEqual([path.name for path in output.iterdir()], ['known.png'])
            SeqIO.write(records[:1], genbank, 'genbank')
            pharokka_plotter.plot(genbank, output)
            self.assertEqual(run.call_count, 1)
            self.assertEqual(list(output.iterdir()), [])

    def test_pharokka_plot_error_is_not_treated_as_empty_annotation(self):
        record = SeqRecord(Seq('ATG' * 10), id='known', name='known', description='')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(FeatureLocation(0, 9), type='CDS',
                                      qualifiers={'product': ['tail protein']})]
        genbank = self.root / 'known.gbk'
        SeqIO.write([record], genbank, 'genbank')
        with patch.object(pharokka_plotter.subprocess, 'run',
                          side_effect=subprocess.CalledProcessError(1, 'pharokka')):
            with self.assertRaises(subprocess.CalledProcessError):
                pharokka_plotter.plot(genbank, self.root / 'failed_plots')

    def test_empty_fasta_and_annotations(self):
        self.fasta.write_text('')
        self.genes.write_text('')
        converter.convert(self.genes, self.fasta, self.gff)
        self.assertEqual(dict(plotter.read_features(self.gff, {})), {})
        args = ['plotter', '--fasta', str(self.fasta), '--genotate-gff', str(self.gff), '--pharokka-gff', str(self.gff), '--output-dir', str(self.root / 'empty_plots')]
        with patch.object(sys, 'argv', args): plotter.main()
        self.assertTrue((self.root / 'empty_plots').is_dir())


if __name__ == '__main__':
    unittest.main()
