import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/prepare_lovis4u_gff.py"
spec = importlib.util.spec_from_file_location("prepare_lovis4u_gff", SCRIPT)
splitter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(splitter)


class PrepareLovis4uGffTests(unittest.TestCase):
    def test_multicontig_split_and_stale_input_removal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "sample" / "pharokka" / "pharokka.gff"
            source.parent.mkdir(parents=True)
            source.write_text(
                "##gff-version 3\n"
                "##sequence-region first 1 9\n"
                "##sequence-region second 1 9\n"
                "first\tPharokka\tCDS\t1\t9\t.\t+\t0\tID=cds1\n"
                "second\tPharokka\tCDS\t1\t9\t.\t+\t0\tID=cds1\n"
                "##FASTA\n>first description\nATGAAATAA\n"
                ">second description\nATGCCCTAA\n"
            )
            output = root / "gff_input"
            output.mkdir()
            (output / "obsolete.gff").write_text("old")

            splitter.prepare([source], output)

            files = sorted(output.glob("*.gff"))
            self.assertEqual(len(files), 2)
            self.assertFalse((output / "obsolete.gff").exists())
            for file in files:
                text = file.read_text()
                self.assertEqual(text.count("##sequence-region"), 1)
                self.assertEqual(text.count("\tCDS\t"), 1)
                self.assertEqual(text.count("\n>"), 1)
                self.assertEqual(text.count("##FASTA"), 1)
                self.assertIn("sample__", file.name)

    def test_missing_sequence_keeps_previous_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "sample" / "pharokka" / "pharokka.gff"
            source.parent.mkdir(parents=True)
            source.write_text(
                "##gff-version 3\n"
                "missing\tPharokka\tCDS\t1\t9\t.\t+\t0\tID=cds1\n"
                "##FASTA\n>other\nATGAAATAA\n"
            )
            output = root / "gff_input"
            output.mkdir()
            (output / "existing.gff").write_text("existing")

            with self.assertRaisesRegex(ValueError, "missing FASTA sequence"):
                splitter.prepare([source], output)
            self.assertEqual((output / "existing.gff").read_text(), "existing")


if __name__ == "__main__":
    unittest.main()
