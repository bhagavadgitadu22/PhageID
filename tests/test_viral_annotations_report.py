"""Checks for the Pharokka viral contig annotation report."""

import ast
import csv
import importlib.util
import re
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "viral_annotations_report", ROOT / "workflow/scripts/viral_annotations_report.py"
)
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def concatenate_reports(reports, output):
    source = (ROOT / "workflow/rules/10_annotation_prodigalgv.smk").read_text()
    section = source.split("rule viral_annotations_report_global:", 1)[1]
    header = ast.literal_eval(re.search(r'header=("[^"]*")', section).group(1))
    if header != "\t".join(report.FIELDS):
        raise AssertionError("Rule header differs from per-sample report header")
    command = ast.literal_eval(re.search(r'shell:\s*(""".*?""")', section, re.S).group(1))
    command = command.replace(":q}", "}").format(
        input=SimpleNamespace(reports=shlex.join([str(path) for path in reports])),
        output=shlex.quote(str(output)),
        params=SimpleNamespace(header=shlex.quote(header)),
        log=shlex.quote(str(output.parent / "global.log")),
    )
    subprocess.run(["bash", "-euo", "pipefail", "-c", command], check=True)


class ViralAnnotationsReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / "viral.fasta"
        self.fasta.write_text(">sample_one\n" + "A" * 100 + "\n>sample_two\n" + "C" * 50 + "\n")
        self.gff = self.root / "pharokka.gff"
        self.gff.write_text(
            "##gff-version 3\n"
            "sample_one\tPharokka\tCDS\t1\t30\t.\t+\t0\tID=cds1\n"
            "sample_one\tPharokka\tCDS\t21\t50\t.\t-\t0\tID=cds2\n"
            "sample_two\tPharokka\ttRNA\t1\t10\t.\t+\t.\tID=trna\n"
            "##FASTA\n>sample_one\nAAAA\n"
        )
        self.functions = self.root / "pharokka_cds_functions.tsv"
        self.functions.write_text(
            "Description\tCount\tcontig\n"
            "CDS\t2\tsample_one\n"
            "tail\t1\tsample_one\n"
            "unknown function\t1\tsample_one\n"
            "CDS\t0\tsample_two\n"
            "tRNAs\t1\tsample_two\n"
        )
        self.output = self.root / "sample.tsv"

    def test_sample_report_counts_unique_cds_bases_and_pivots_categories(self):
        report.sample_report("sample", self.fasta, self.gff, self.functions, self.output)
        result = rows(self.output)
        self.assertEqual(list(result[0]), report.FIELDS)
        self.assertEqual([row["Viral_contig"] for row in result], ["sample_one", "sample_two"])
        self.assertEqual(result[0]["Coding_density_percent"], "50.00")
        self.assertEqual(result[0]["CDS"], "2")
        self.assertEqual(result[0]["tail"], "1")
        self.assertEqual(result[0]["lysis"], "0")
        self.assertEqual(result[1]["Coding_density_percent"], "0.00")
        self.assertEqual(result[1]["tRNAs"], "1")

    def test_global_rule_concatenates_sample_rows_once(self):
        report.sample_report("sample", self.fasta, self.gff, self.functions, self.output)
        second_report = self.root / "second.tsv"
        second = rows(self.output)[0]
        second["Viral_contig"] = "another_contig"
        second["Sample"] = "another"
        report.write_tsv(second_report, [second])
        global_output = self.root / "global.tsv"
        concatenate_reports([self.output, second_report], global_output)
        self.assertEqual([row["Viral_contig"] for row in rows(global_output)],
                         ["sample_one", "sample_two", "another_contig"])
        self.assertEqual(global_output.read_text().count("Viral_contig\tSample"), 1)

    def test_global_rule_without_samples_has_header(self):
        output = self.root / "global.tsv"
        concatenate_reports([], output)
        self.assertEqual(output.read_text().strip().split("\t"), report.FIELDS)

    def test_rejects_unknown_function_contig(self):
        self.functions.write_text("Description\tCount\tcontig\nCDS\t1\tother\n")
        with self.assertRaisesRegex(ValueError, "unknown contig"):
            report.sample_report("sample", self.fasta, self.gff, self.functions, self.output)


if __name__ == "__main__":
    unittest.main()
