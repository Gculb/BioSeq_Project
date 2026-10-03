import json
import os
import tempfile
import unittest

from bioseq.pipeline import run_pipeline
from bioseq.qc import analyze_fastq_quality
from bioseq.validation import detect_sequence_format, validate_sequence_file


class SequenceValidationTests(unittest.TestCase):
    def test_valid_fasta(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", encoding="utf-8", delete=False
        ) as file_handle:
            file_handle.write(">seq1 description\nACGTN\n>seq2\nGGCC\n")
            filename = file_handle.name

        try:
            result = validate_sequence_file(filename)
        finally:
            os.unlink(filename)

        self.assertTrue(result["valid"])
        self.assertEqual(result["format"], "fasta")
        self.assertEqual(result["record_count"], 2)

    def test_invalid_fasta_reports_empty_sequence_and_invalid_bases(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fa", encoding="utf-8", delete=False
        ) as file_handle:
            file_handle.write(">empty\n>bad\nAC-Z\n")
            filename = file_handle.name

        try:
            result = validate_sequence_file(filename)
        finally:
            os.unlink(filename)

        self.assertFalse(result["valid"])
        self.assertEqual(len(result["errors"]), 2)

    def test_valid_fastq_and_invalid_quality_range(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fastq", encoding="utf-8", delete=False
        ) as file_handle:
            file_handle.write("@read\nACGT\n+\nII I\n")
            filename = file_handle.name

        try:
            result = validate_sequence_file(filename)
        finally:
            os.unlink(filename)

        self.assertFalse(result["valid"])
        self.assertTrue(any("Phred+33" in error for error in result["errors"]))

    def test_unsupported_extension_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            detect_sequence_format("sequences.txt")


class FastqQualityControlTests(unittest.TestCase):
    def test_quality_metrics_are_calculated(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fastq", encoding="utf-8", delete=False
        ) as file_handle:
            file_handle.write(
                "@high\nACGTN\n+\nIIIII\n"
                "@low\nGGCC\n+\n!!!!\n"
            )
            filename = file_handle.name

        try:
            metrics = analyze_fastq_quality(filename)
        finally:
            os.unlink(filename)

        self.assertEqual(metrics["read_count"], 2)
        self.assertEqual(metrics["total_bases"], 9)
        self.assertEqual(metrics["read_length_distribution"], {4: 1, 5: 1})
        self.assertAlmostEqual(metrics["mean_quality"], 200 / 9)
        self.assertAlmostEqual(metrics["gc_content"], 600 / 9)
        self.assertEqual(metrics["ambiguous_bases"], 1)
        self.assertEqual(metrics["low_quality_reads"], 1)
        self.assertEqual(metrics["low_quality_percentage"], 50)


class PipelineTests(unittest.TestCase):
    def test_pipeline_validates_analyzes_and_writes_both_reports(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "reads.fastq")
            results_dir = os.path.join(temp_dir, "results")
            with open(input_file, "w", encoding="utf-8") as file_handle:
                file_handle.write("@read\nACGT\n+\nIIII\n")

            report = run_pipeline(input_file, results_dir=results_dir)

            self.assertTrue(report["validation"]["valid"])
            self.assertEqual(report["format"], "fastq")
            self.assertEqual(report["metrics"]["read_count"], 1)
            for report_path in report["report_files"].values():
                self.assertTrue(os.path.isfile(report_path))

            with open(report["report_files"]["json"], encoding="utf-8") as file_handle:
                saved_report = json.load(file_handle)
            self.assertEqual(saved_report["metrics"]["read_count"], 1)

    def test_pipeline_stops_on_invalid_input(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "bad.fasta")
            with open(input_file, "w", encoding="utf-8") as file_handle:
                file_handle.write(">sequence\nAC?T\n")

            with self.assertRaisesRegex(ValueError, "Input validation failed"):
                run_pipeline(input_file, results_dir=temp_dir)


if __name__ == "__main__":
    unittest.main()
