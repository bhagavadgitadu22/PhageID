import importlib.util
import tempfile
import unittest
from pathlib import Path
from Bio import SeqIO

spec = importlib.util.spec_from_file_location('dereplication_mapping', Path(__file__).resolve().parents[1] / 'workflow/scripts/dereplication_mapping.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class MappingTests(unittest.TestCase):
    def test_membership_and_single_representative_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p/'reps.fa').write_text('>other_contig\nACGT\n>s_one_contig2\nAAAA\n')
            (p/'s.fa').write_text('>s_one_contig1\nACGT\n>s_one_contig2\nAAAA\n>s_one_contig3\nACGT\n')
            (p/'other.fa').write_text('>other_contig\nACGT\n')
            (p/'clusters.tsv').write_text('other_contig\tother_contig,s_one_contig1,s_one_contig3\ns_one_contig2\ts_one_contig2\n')
            module.prepare(p/'reps.fa', p/'clusters.tsv', [('s_one', p/'s.fa'), ('other', p/'other.fa')], p/'out')
            self.assertEqual(module.read_pairs(p/'out/sample_representatives.tsv'), [('other', 'other_contig'), ('s_one', 'other_contig'), ('s_one', 's_one_contig2')])
            for representative in ['other_contig', 's_one_contig2']:
                records = list(SeqIO.parse(p/'out/references'/f'{representative}.fna', 'fasta'))
                self.assertEqual([r.id for r in records], [representative])
            (p/'clusters.tsv').write_text('other_contig\tother_contig,unknown\n')
            with self.assertRaisesRegex(ValueError, 'Unknown or duplicated cluster member'):
                module.prepare(p/'reps.fa', p/'clusters.tsv', [('s_one', p/'s.fa'), ('other', p/'other.fa')], p/'invalid')

    def test_empty_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            for name in ['reps.fa', 'clusters.tsv', 'sample.fa']: (p/name).write_text('')
            module.prepare(p/'reps.fa', p/'clusters.tsv', [('sample', p/'sample.fa')], p/'out')
            self.assertEqual(module.read_pairs(p/'out/sample_representatives.tsv'), [])


if __name__ == '__main__': unittest.main()
