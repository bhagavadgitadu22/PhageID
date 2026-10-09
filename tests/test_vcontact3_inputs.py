import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prepare_vcontact3', ROOT/'workflow/scripts/prepare_vcontact3.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class VcontactInputsTests(unittest.TestCase):
    def test_only_representatives_and_matching_unique_protein_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/'representatives.fa').write_text('>sample_a_contig\nATGAAATAA\n>sample_b_contig\nATGAAATAA\n')
            records=[]
            for name in ['sample_a_contig','discarded','sample_b_contig']:
                record=SeqRecord(Seq('ATGAAATAA'),id=name,name=name,description='')
                record.annotations['molecule_type']='DNA'
                record.features=[SeqFeature(FeatureLocation(0,9),type='CDS',qualifiers={'locus_tag':['same_tag'],'translation':['MK']})]
                records.append(record)
            SeqIO.write(records,p/'annotations.gbk','genbank')
            module.prepare(p/'representatives.fa',[p/'annotations.gbk'],p/'proteins.faa',p/'map.tsv',p/'lengths.tsv')
            proteins=list(SeqIO.parse(p/'proteins.faa','fasta'))
            with open(p/'map.tsv') as handle:
                mappings=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r.id for r in proteins],[r['protein_id'] for r in mappings])
            self.assertEqual([r['genome_id'] for r in mappings],['sample_a_contig','sample_b_contig'])
            self.assertEqual(len({r.id for r in proteins}),2)
            self.assertEqual([str(r.seq) for r in proteins],['MK','MK'])
            with open(p/'lengths.tsv') as handle:
                lengths=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['length'] for r in lengths],['9','9'])
            SeqIO.write(records[:1],p/'missing.gbk','genbank')
            with self.assertRaisesRegex(ValueError,'Missing representative annotations'):
                module.prepare(p/'representatives.fa',[p/'missing.gbk'],p/'bad.faa',p/'bad.tsv',p/'bad_lengths.tsv')

    def test_empty_representatives(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/'empty.fa').write_text('')
            module.prepare(p/'empty.fa',[],p/'proteins.faa',p/'map.tsv',p/'lengths.tsv')
            self.assertEqual((p/'proteins.faa').read_text(),'')
            self.assertEqual((p/'map.tsv').read_text(),'protein_id\tgenome_id\n')
