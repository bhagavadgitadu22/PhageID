import importlib.util
import unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

spec = importlib.util.spec_from_file_location('correction', Path(__file__).resolve().parents[1] / 'workflow/scripts/correct_phage_contigs.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def hit(qstart, qend, sstart, send, identity=100):
    return dict(qstart=qstart, qend=qend, sstart=sstart, send=send,
                length=abs(qend-qstart)+1, pident=identity)


class ConcatemerTests(unittest.TestCase):
    unit = 'ACGTTGCAAG' * 10

    def correct(self, sequence, hits):
        return module.correct_concatemer(SeqRecord(Seq(sequence), id='contig'), hits, 95, 90, 90)

    def test_two_three_and_four_copies(self):
        for copies in (2, 3, 4):
            with self.subTest(copies=copies):
                size = len(self.unit) * copies
                record, report = self.correct(self.unit * copies, [hit(1,size-100,101,size)])
                self.assertEqual(str(record.seq), self.unit)
                self.assertEqual(report['repeat_unit_size'], 100)
                self.assertEqual(report['num_copies'], copies)
                self.assertEqual(report['status'], 'corrected')

    def test_imperfect_copy_within_identity_threshold(self):
        imperfect = 'T' + self.unit[1:]
        record, report = self.correct(self.unit + imperfect + self.unit, [hit(1,200,101,300,99)])
        self.assertEqual(str(record.seq), self.unit)
        self.assertEqual(report['num_copies'], 3)

    def test_whole_sequence_must_support_period(self):
        sequence = self.unit + 'T' * 100 + self.unit
        record, report = self.correct(sequence, [hit(1,200,101,300)])
        self.assertEqual(str(record.seq), sequence)
        self.assertEqual(report['status'], 'repeats_unclear')

    def test_internal_and_terminal_repeats_remain_unchanged(self):
        cases = [
            (self.unit + 'T'*100 + self.unit, [hit(1,100,201,300)]),
            ('T'*100 + self.unit*2 + 'G'*100, [hit(101,200,201,300)]),
        ]
        for sequence, hits in cases:
            with self.subTest(sequence=sequence):
                record, report = self.correct(sequence, hits)
                self.assertEqual(str(record.seq), sequence)
                self.assertFalse(report['concatemer_detected'])

    def test_partial_copy_is_not_trimmed(self):
        sequence = self.unit*2 + self.unit[:50]
        record, report = self.correct(sequence, [hit(1,150,101,250)])
        self.assertEqual(str(record.seq), sequence)
        self.assertEqual(report['status'], 'repeats_unclear')

    def test_inverted_and_no_repeats(self):
        for hits in ([], [hit(1,200,300,101)]):
            record, report = self.correct(self.unit*3, hits)
            self.assertEqual(len(record.seq), 300)
            self.assertFalse(report['concatemer_detected'])

    def test_duplicate_hits_do_not_inflate_coverage(self):
        sequence = self.unit*3
        record, report = self.correct(sequence, [hit(1,30,101,130)]*10)
        self.assertEqual(str(record.seq), sequence)
        self.assertEqual(report['unique_repeat_coverage_bp'],60)


if __name__ == '__main__': unittest.main()
