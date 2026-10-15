import importlib.util
import tempfile
import unittest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prepare_vcontact3', ROOT/'workflow/scripts/prepare_vcontact3.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class VcontactInputsTests(unittest.TestCase):
    @staticmethod
    def write_genbank(path, identifier, sequence='ATGAAATAA'):
        path.parent.mkdir(parents=True, exist_ok=True)
        record = SeqRecord(Seq(sequence), id=identifier, name=identifier, description='')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(FeatureLocation(0, len(sequence)), type='CDS',
                                      qualifiers={'translation': ['MK']})]
        SeqIO.write([record], path, 'genbank')

    def test_only_representatives_and_matching_unique_protein_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/'representatives.fa').write_text('>sample_a_contig\nATGAAATAA\n>sample_b_contig\nATGAAATAA\n')
            records=[]
            for name in ['sample_a_contig','discarded','sample_b_contig']:
                record=SeqRecord(Seq('ATGAAATAA'),id=name,name=name,description='')
                record.annotations['molecule_type']='DNA'
                qualifiers={'locus_tag':['same_tag'],'translation':['MK']}
                if name == 'sample_a_contig':
                    qualifiers['product']=['capsid protein']
                record.features=[SeqFeature(FeatureLocation(0,9),type='CDS',qualifiers=qualifiers)]
                records.append(record)
            SeqIO.write(records,p/'annotations.gbk','genbank')
            module.prepare(p/'representatives.fa',[p/'annotations.gbk'],p/'proteins.faa',p/'map.parquet',p/'lengths.parquet')
            proteins=list(SeqIO.parse(p/'proteins.faa','fasta'))
            mappings=pd.read_parquet(p/'map.parquet')
            self.assertEqual(mappings.columns.tolist(),['protein_id','keywords','genome_id'])
            self.assertEqual([r.id for r in proteins],mappings['protein_id'].tolist())
            self.assertEqual(mappings['genome_id'].tolist(),['sample_a_contig','sample_b_contig'])
            self.assertEqual(mappings['keywords'].tolist(),['capsid protein','None'])
            overview=mappings.replace('None',pd.NA).groupby('genome_id').nunique()
            overview['genome_id']=overview.index.tolist()
            self.assertEqual(len(overview.columns),3)
            self.assertEqual(len({r.id for r in proteins}),2)
            self.assertEqual([str(r.seq) for r in proteins],['MK','MK'])
            lengths=pd.read_parquet(p/'lengths.parquet')
            self.assertEqual(lengths['genome_id'].tolist(),['sample_a_contig','sample_b_contig'])
            self.assertEqual(lengths['length'].tolist(),[9,9])
            SeqIO.write(records[:1],p/'missing.gbk','genbank')
            with self.assertRaisesRegex(ValueError,'Missing representative annotations'):
                module.prepare(p/'representatives.fa',[p/'missing.gbk'],p/'bad.faa',p/'bad.parquet',p/'bad_lengths.parquet')

    def test_empty_representatives(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/'empty.fa').write_text('')
            module.prepare(p/'empty.fa',[],p/'proteins.faa',p/'map.parquet',p/'lengths.parquet')
            self.assertEqual((p/'proteins.faa').read_text(),'')
            self.assertEqual(pd.read_parquet(p/'map.parquet').columns.tolist(),['protein_id','keywords','genome_id'])
            self.assertEqual(pd.read_parquet(p/'lengths.parquet').columns.tolist(),['genome_id','length'])

    def test_directory_and_quoted_glob_find_same_files_once(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/'representatives.fa').write_text('>first\nATGAAATAA\n>second\nATGAAATAA\n')
            self.write_genbank(p/'sample_a'/'pharokka'/'pharokka.gbk', 'first')
            self.write_genbank(p/'sample_b'/'pharokka'/'pharokka.gbk', 'second')
            pattern = str(p/'*'/'pharokka'/'*.gbk')
            module.prepare(p/'representatives.fa', [pattern, p/'sample_a'],
                           p/'proteins.faa', p/'map.parquet', p/'lengths.parquet')
            self.assertEqual(pd.read_parquet(p/'map.parquet')['genome_id'].tolist(),
                             ['first', 'second'])

    def test_unmatched_pattern_and_missing_id_report_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/'representatives.fa').write_text('>first\nATGAAATAA\n>second\nATGAAATAA\n')
            annotation = p/'pharokka.gbk'
            self.write_genbank(annotation, 'first')
            with self.assertRaisesRegex(ValueError, 'No GenBank files match'):
                module.prepare(p/'representatives.fa', [str(p/'absent'/'*.gbk')],
                               p/'proteins.faa', p/'map.parquet', p/'lengths.parquet')
            with self.assertRaisesRegex(ValueError, 'Missing representative annotations.*second.*searched 1 GenBank files') as error:
                module.prepare(p/'representatives.fa', [annotation],
                               p/'proteins.faa', p/'map.parquet', p/'lengths.parquet')
            self.assertIn(str(annotation), str(error.exception))
            self.assertFalse((p/'proteins.faa').exists())

    def test_same_id_with_different_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory)
            (p/'representatives.fa').write_text('>first\nATGAAATAA\n')
            annotation = p/'pharokka.gbk'
            self.write_genbank(annotation, 'first', 'ATGCCCTAA')
            with self.assertRaisesRegex(ValueError, 'Annotation sequence differs for first'):
                module.prepare(p/'representatives.fa', [annotation],
                               p/'proteins.faa', p/'map.parquet', p/'lengths.parquet')
