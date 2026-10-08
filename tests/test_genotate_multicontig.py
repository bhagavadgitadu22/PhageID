"""Regression checks for multi-contig conversion and plot routing."""
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'workflow' / 'scripts' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


converter = load('reformat_pharokka_genotate')
plotter = load('plotting_pharokka_genotate')


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
        self.assertEqual(calls, [{'a': 100}, {'b': 200}, {'no_genes': 4}])
        self.assertEqual(len(list(output.glob('*.png'))), 3)
        self.assertEqual(len(arrows), 4)  # Each CDS once in each annotation source.

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
