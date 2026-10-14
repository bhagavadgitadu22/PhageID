import importlib.util
import tempfile
import unittest
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

spec = importlib.util.spec_from_file_location(
    "genbank_to_gff", Path(__file__).resolve().parents[1] / "workflow/scripts/genbank_to_gff.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ConversionTests(unittest.TestCase):
    def test_qualifiers_coordinates_and_joined_reverse_cds(self):
        record = SeqRecord(Seq("A" * 30), id="contig", annotations={"molecule_type": "DNA"})
        record.features = [SeqFeature(
            CompoundLocation([SimpleLocation(20, 30, strand=-1), SimpleLocation(0, 5, strand=-1)]),
            type="CDS", qualifiers={"locus_tag": ["gene_1"], "codon_start": ["2"],
                                    "product": ["capsid; protein, A=B"], "note": ["one", "two"]})]
        other = SeqRecord(Seq("AAA"), id="other", annotations={"molecule_type": "DNA"})
        other.features = [SeqFeature(SimpleLocation(0, 3, strand=1), type="CDS",
                                     qualifiers={"locus_tag": ["gene_2"]})]
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / "phold.gbk", Path(tmp) / "phold.gff"
            SeqIO.write([record, other], source, "genbank")
            module.convert(source, target)
            text = target.read_text()
        self.assertIn("##sequence-region contig 1 30", text)
        self.assertIn("##sequence-region other 1 3", text)
        rows = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
        self.assertEqual([row[3:8] for row in rows],
                         [["21", "30", ".", "-", "1"], ["1", "5", ".", "-", "0"],
                          ["1", "3", ".", "+", "0"]])
        attributes = dict(field.split("=", 1) for field in rows[0][8].split(";"))
        self.assertEqual(attributes["locus_tag"], "gene_1")
        self.assertEqual(unquote(attributes["product"]), "capsid; protein, A=B")
        self.assertEqual(attributes["note"], "one,two")
        self.assertEqual(rows[0][8], rows[1][8])
        self.assertNotEqual(rows[0][8].split("ID=")[1], rows[2][8].split("ID=")[1])


if __name__ == "__main__":
    unittest.main()
