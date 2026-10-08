import importlib.util
import tempfile
import unittest
from pathlib import Path
from Bio import SeqIO

spec = importlib.util.spec_from_file_location('filter_annotations', Path(__file__).resolve().parents[1] / 'workflow/scripts/filter_dereplicated_annotations.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FilterTests(unittest.TestCase):
    def test_only_representative_features_and_sequences_are_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            reps = p/'reps.fa'
            reps.write_text('>sample_b_rep\nACGTACGT\n>gene_less\nAAA\n')
            first = p/'a.gff'
            first.write_text('##gff-version 3\n##sequence-region sample_a_member 1 8\nsample_a_member\tphold\tCDS\t1\t6\t.\t+\t0\tID=removed\n##FASTA\n>sample_a_member\nACGTACGT\n')
            second = p/'b.gff'
            feature = 'sample_b_rep\tphold\tCDS\t1\t6\t.\t-\t0\tID=kept;empathi_Annotation=lysis;product=capsid%20protein\n'
            second.write_text('##gff-version 3\n##sequence-region sample_b_rep 1 8\n'+feature+'##FASTA\n>sample_b_rep\nTTTTTTTT\n')
            out = p/'out.gff'
            module.filter_annotations([first,second], reps, out)
            text = out.read_text()
            self.assertIn(feature, text)
            self.assertNotIn('sample_a_member', text)
            self.assertEqual(text.count('##gff-version'),1)
            self.assertEqual(text.count('##FASTA'),1)
            with out.open() as handle:
                for line in handle:
                    if line.startswith('##FASTA'): break
                self.assertEqual([(r.id,str(r.seq)) for r in SeqIO.parse(handle,'fasta')], [('sample_b_rep','ACGTACGT'),('gene_less','AAA')])

    def test_merge_all_samples_with_embedded_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            annotations, fastas = [], []
            for contig, sequence in [('a', 'ACGTACGT'), ('b', 'TTTTAAAA')]:
                fasta = p / f'{contig}.fa'
                fasta.write_text(f'>{contig}\n{sequence}\n')
                gff = p / f'{contig}.gff'
                gff.write_text(f'##gff-version 3\n{contig}\tphold\tCDS\t1\t6\t.\t+\t0\tID={contig}_CDS;empathi_Annotation=lysis\n##FASTA\n>{contig}\n{sequence}\n')
                annotations.append(gff)
                fastas.append(fasta)
            module.filter_annotations(annotations, fastas, p/'merged.gff')
            text = (p/'merged.gff').read_text()
            features, sequences = text.split('##FASTA\n')
            self.assertIn('a\tphold\tCDS', features)
            self.assertIn('b\tphold\tCDS', features)
            self.assertEqual(text.count('##gff-version'), 1)
            self.assertIn('>a\nACGTACGT', sequences)
            self.assertIn('>b\nTTTTAAAA', sequences)

    def test_empty_representatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); (p/'reps.fa').write_text('')
            module.filter_annotations([],p/'reps.fa',p/'out.gff')
            self.assertEqual((p/'out.gff').read_text(),'##gff-version 3\n')

    def test_invalid_coordinates_do_not_replace_existing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); (p/'reps.fa').write_text('>rep\nACGT\n')
            (p/'in.gff').write_text('rep\tx\tCDS\t1\t9\t.\t+\t0\tID=x\n')
            (p/'out.gff').write_text('existing')
            with self.assertRaisesRegex(ValueError,'coordinates outside'):
                module.filter_annotations([p/'in.gff'],p/'reps.fa',p/'out.gff')
            self.assertEqual((p/'out.gff').read_text(),'existing')


if __name__ == '__main__': unittest.main()
