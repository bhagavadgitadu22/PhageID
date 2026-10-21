"""Representative ranking and directional ANI behavior for viral dereplication."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/ani_clust.py"
QUALITY_FIELDS = ["Viral_contig", "CheckV_quality", "CheckV_completeness",
                  "CheckV_contamination"]


class AniClusterTests(unittest.TestCase):
    def run_cluster(self, contigs, quality, alignments):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with open(root / "contigs.fna", "w") as handle:
                for name, length in contigs.items():
                    handle.write(f">{name}\n" + "A" * length + "\n")
            with open(root / "quality.tsv", "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(QUALITY_FIELDS)
                for name, values in quality.items():
                    writer.writerow([name, *values])
            with open(root / "ani.tsv", "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["qname", "tname", "num_alns", "ani", "qcov", "tcov"])
                for query, target, qcov, tcov in alignments:
                    writer.writerow([query, target, 1, 99, qcov, tcov])
            result = subprocess.run(
                [sys.executable, "-S", str(SCRIPT), "--fna", str(root / "contigs.fna"),
                 "--ani", str(root / "ani.tsv"), "--quality-report", str(root / "quality.tsv"),
                 "--out", str(root / "clusters.tsv"), "--min_ani", "95",
                 "--min_tcov", "85", "--min_qcov", "0"],
                capture_output=True, text=True)
            clusters = (root / "clusters.tsv").read_text().splitlines() if result.returncode == 0 else []
            return result, clusters

    def test_quality_category_beats_length_and_numeric_score(self):
        result, clusters = self.run_cluster(
            {"medium": 120, "high": 100},
            {"medium": ("Medium-quality", "95", "0"),
             "high": ("High-quality", "60", "10")},
            [("high", "medium", 100, 85), ("medium", "high", 85, 100)])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(clusters, ["high\thigh,medium"])

    def test_score_then_missing_values_then_length_and_id(self):
        result, clusters = self.run_cluster(
            {"contaminated": 120, "clean": 100, "unknown": 200},
            {"contaminated": ("High-quality", "95", "5"),
             "clean": ("High-quality", "90", "0"),
             "unknown": ("High-quality", "NA", "0")},
            [("clean", "contaminated", 100, 100),
             ("clean", "unknown", 100, 100)])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(clusters, ["clean\tclean,contaminated,unknown"])

        result, clusters = self.run_cluster(
            {"zeta": 100, "alpha": 100},
            {"zeta": ("Not-determined", "NA", "NA"),
             "alpha": ("Not-determined", "NA", "NA")},
            [("alpha", "zeta", 100, 100)])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(clusters, ["alpha\talpha,zeta"])

    def test_quality_first_keeps_directional_member_coverage_threshold(self):
        result, clusters = self.run_cluster(
            {"medium": 150, "high": 100},
            {"medium": ("Medium-quality", "70", "0"),
             "high": ("High-quality", "90", "0")},
            [("high", "medium", 100, 70), ("medium", "high", 70, 100)])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(clusters, ["high\thigh", "medium\tmedium"])

    def test_missing_quality_row_fails_instead_of_silent_length_fallback(self):
        result, _ = self.run_cluster(
            {"one": 100}, {}, [])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Missing CheckV quality for contigs: one", result.stderr)


if __name__ == "__main__":
    unittest.main()
