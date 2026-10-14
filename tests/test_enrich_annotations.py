import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "enrich_annotations", Path(__file__).resolve().parents[1] / "workflow/scripts/enrich_annotations.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class EnrichmentTests(unittest.TestCase):
    def test_empty_results_preserve_input_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, csv, output = [Path(tmp) / name for name in ("input.gff", "hits.csv", "out.gff")]
            source.write_bytes(b"##gff-version 3\ncontig\tphold\tCDS\t1\t9\t.\t+\t0\tlocus_tag=gene_1\n")
            for contents in ("", "feature_type,locus_tag,hit\n", "feature_type,locus_tag,hit\n\n,,\n"):
                with self.subTest(contents=contents), patch.object(module.subprocess, "run") as run:
                    csv.write_text(contents)
                    module.enrich(source, csv, output, "antidefensefinder_")
                    run.assert_not_called()
                    self.assertEqual(output.read_bytes(), source.read_bytes())

    def test_nonempty_results_do_not_hide_tool_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv = Path(tmp) / "hits.csv"
            csv.write_text('feature_type,locus_tag,hit\nCDS,gene_1,"value, with comma"\n')
            with patch.object(module.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "thebigbam")) as run:
                with self.assertRaises(subprocess.CalledProcessError):
                    module.enrich("input.gff", csv, Path(tmp) / "out.gff", "empathi_")
                self.assertTrue(run.call_args.kwargs["check"])


if __name__ == "__main__":
    unittest.main()
