"""Manual FASTA overrides preserve correction provenance and validate IDs."""

import csv
import importlib.util
import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/finalize_viral_assembly.py"
SPEC = importlib.util.spec_from_file_location("finalize_viral_assembly", SCRIPT)
finalize = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(finalize)


class ManualViralCorrectionsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.auto = self.root / "automatic.fasta"
        self.auto.write_text(">sample_a\nAAAACCCC\n>sample_b\nGGGG\n")
        self.concatemer = self.root / "concatemer.csv"
        self.concatemer.write_text(
            "contig_id,original_length,corrected_length\n"
            "sample_a,12,10\nsample_b,4,4\n")
        self.dtr = self.root / "dtr.csv"
        self.dtr.write_text(
            "contig_id,original_length,corrected_length\n"
            "sample_a,10,8\nsample_b,4,4\n")
        self.manual = self.root / "manual.fasta"
        self.output = self.root / "published/final.fasta"
        self.report = self.root / "published/manual.tsv"

    def publish(self, manual=None):
        finalize.publish(self.auto, self.concatemer, self.dtr, manual, self.output,
                         self.root / "published/concatemer.csv",
                         self.root / "published/dtr.csv", self.report)

    def report_rows(self):
        with self.report.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_override_changes_one_contig_and_preserves_automatic_dtr_lengths(self):
        self.manual.write_text(">sample_a\nAAAA\n>sample_b\nGGGG\n")
        self.publish(self.manual)
        self.assertEqual(self.output.read_text(), self.manual.read_text())
        rows = self.report_rows()
        self.assertEqual((rows[0]["Auto_corrected_bp"], rows[0]["Final_bp"],
                          rows[0]["Manual_correction"], rows[0]["Manual_length_change_bp"]),
                         ("8", "4", "Yes", "-4"))
        self.assertEqual(rows[1]["Manual_correction"], "No")
        self.assertEqual((self.root / "published/dtr.csv").read_text(), self.dtr.read_text())

    def test_without_override_records_no_changes(self):
        self.publish()
        self.assertEqual(self.output.read_text(), self.auto.read_text())
        self.assertEqual([row["Manual_correction"] for row in self.report_rows()], ["No", "No"])

    def test_missing_or_reordered_contigs_do_not_replace_published_output(self):
        self.output.parent.mkdir()
        self.output.write_text("previous result")
        for contents in (">sample_a\nAAAA\n", ">sample_b\nGGGG\n>sample_a\nAAAA\n"):
            self.manual.write_text(contents)
            with self.assertRaisesRegex(ValueError, "exactly the automatic contig IDs"):
                self.publish(self.manual)
            self.assertEqual(self.output.read_text(), "previous result")

    def test_empty_manual_contig_identifies_the_record(self):
        self.manual.write_text(">sample_a\nAAAA\n>sample_b\n")
        with self.assertRaisesRegex(ValueError, "Empty sequence for contig sample_b"):
            self.publish(self.manual)
        self.assertFalse(self.output.exists())

    def test_rule_copies_or_recomputes_checkv_according_to_final_sequence(self):
        source = (SCRIPT.parent.parent / "rules/04_assembly_selection.smk").read_text()
        command = source.split("rule finalize_viral_assembly:", 1)[1].split('"""', 1)[1].split('"""', 1)[0]
        quality = self.root / "automatic_quality.tsv"
        quality.write_text("automatic quality\n")
        final_quality = self.root / "published/checkv/quality_summary.tsv"
        mock_bin = self.root / "bin"
        mock_bin.mkdir()
        mock_checkv = mock_bin / "checkv"
        mock_checkv.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\nfrom pathlib import Path\n"
            "(Path(sys.argv[-1]) / 'quality_summary.tsv').write_text('recomputed quality\\n')\n")
        mock_checkv.chmod(0o755)

        def run(manual_args):
            quote = lambda value: shlex.quote(str(value))
            rendered = command.replace(":q}", "}").format(
                params=SimpleNamespace(script=quote(SCRIPT), manual_args=shlex.join(manual_args)),
                input=SimpleNamespace(automatic=quote(self.auto), concatemer=quote(self.concatemer),
                                      dtr=quote(self.dtr), quality=quote(quality), db=quote(self.root)),
                output=SimpleNamespace(corrected=quote(self.output),
                                       concatemer_report=quote(self.root / "published/concatemer.csv"),
                                       dtr_report=quote(self.root / "published/dtr.csv"),
                                       manual_report=quote(self.report), checkv_quality=quote(final_quality)),
                threads=1, log=quote(self.root / "finalize.log"))
            environment = {**os.environ, "PATH": str(mock_bin) + os.pathsep + os.environ["PATH"]}
            return subprocess.run(["bash", "-euo", "pipefail", "-c", rendered],
                                  capture_output=True, text=True, env=environment)

        result = run([])
        self.assertEqual(result.returncode, 0, (result.stderr, (self.root / "finalize.log").read_text()))
        self.assertEqual(final_quality.read_text(), "automatic quality\n")
        self.manual.write_text(">sample_a\nAAAA\n>sample_b\nGGGG\n")
        result = run(["--manual", str(self.manual)])
        self.assertEqual(result.returncode, 0, (result.stderr, (self.root / "finalize.log").read_text()))
        self.assertEqual(final_quality.read_text(), "recomputed quality\n")


if __name__ == "__main__":
    unittest.main()
