import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/reformat_feature_annotations.py"
spec = importlib.util.spec_from_file_location("reformat_feature_annotations", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CheckamgReformatTests(unittest.TestCase):
    def run_reformatter(self, directory, raw_headers, processed_headers, result_ids):
        raw = directory / "prodigal-gv_raw.faa"
        processed = directory / "prodigal-gv.faa"
        annotations = directory / "final_results.tsv"
        output = directory / "reformatted.csv"
        raw.write_text("".join(f">{header}\nM\n" for header in raw_headers))
        processed.write_text("".join(f">{header}\nM\n" for header in processed_headers))
        annotations.write_text(
            "protein_id\tclassification\n"
            + "".join(f"{identifier}\tAMG\n" for identifier in result_ids)
        )
        argv = [
            str(SCRIPT), "--input", str(annotations), "--output", str(output),
            "--checkamg-pharokka-fasta", str(processed),
            "--checkamg-raw-fasta", str(raw),
        ]
        with patch("sys.argv", argv):
            module.main()
        with output.open(newline="") as handle:
            return list(csv.DictReader(handle))

    def test_pairs_ids_by_fasta_order_across_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = self.run_reformatter(
                Path(tmp),
                ["sample_contig_1_1", "sample_contig_1_10", "sample_contig_2_1"],
                ["ABC_CDS_0001", "ABC_CDS_0010", "XYZ_CDS_0001"],
                ["sample_contig_2_1", "sample_contig_1_10"],
            )
        self.assertEqual([row["locus_tag"] for row in rows], ["XYZ_CDS_0001", "ABC_CDS_0010"])
        self.assertEqual(
            [row["checkamg_id"] for row in rows],
            ["sample_contig_2_1", "sample_contig_1_10"],
        )

    def test_rejects_unequal_fasta_record_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "record counts differ"):
                self.run_reformatter(
                    Path(tmp), ["sample_contig_1_1", "sample_contig_1_2"],
                    ["ABC_CDS_0001"], ["sample_contig_1_1"],
                )

    def test_rejects_result_id_missing_from_raw_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "no Pharokka protein for CheckAMG ID"):
                self.run_reformatter(
                    Path(tmp), ["sample_contig_1_1"],
                    ["ABC_CDS_0001"], ["sample_contig_1_2"],
                )


if __name__ == "__main__":
    unittest.main()
