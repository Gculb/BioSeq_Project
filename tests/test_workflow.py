import json
import gzip
import os
import tempfile
import unittest
from unittest.mock import patch

from bioseq.pipeline import run_pipeline
from bioseq.qc import analyze_fastq_quality
from bioseq.validation import detect_sequence_format, validate_sequence_file
from bioseq.variants import analyze_vcf


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

    def test_variant_and_alignment_formats_are_detected(self):
        self.assertEqual(detect_sequence_format("calls.vcf"), "vcf")
        self.assertEqual(detect_sequence_format("calls.vcf.gz"), "vcf")
        self.assertEqual(detect_sequence_format("reads.bam"), "bam")
        self.assertEqual(detect_sequence_format("reads.cram"), "cram")


class VariantAnalysisTests(unittest.TestCase):
    def test_analyze_vcf_counts_variants_filters_and_genotypes(self):
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n"
            "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:12\t1/1:8\n"
            "chr1\t20\trs2\tAT\tA,ATT\t.\tLowQual\t.\tGT\t./.\t0/2\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "calls.vcf.gz")
            with gzip.open(filename, "wt", encoding="utf-8") as file_handle:
                file_handle.write(vcf_text)

            result = analyze_vcf(filename)

        self.assertTrue(result["validation"]["valid"])
        self.assertEqual(result["validation"]["record_count"], 2)
        self.assertEqual(result["metrics"]["snps"], 1)
        self.assertEqual(result["metrics"]["indels"], 2)
        self.assertEqual(result["metrics"]["pass_records"], 1)
        self.assertEqual(result["metrics"]["filtered_records"], 1)
        self.assertEqual(result["metrics"]["samples"], ["sample1", "sample2"])
        self.assertEqual(result["metrics"]["called_genotypes"], 3)
        self.assertEqual(result["metrics"]["missing_genotypes"], 1)
        self.assertEqual(result["metrics"]["heterozygous_genotypes"], 2)
        self.assertEqual(result["metrics"]["homozygous_alternate_genotypes"], 1)
        self.assertEqual(result["metrics"]["contig_count"], 1)
        self.assertEqual(result["metrics"]["mean_qual"], 50)

    def test_invalid_vcf_fails_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "calls.vcf")
            with open(filename, "w", encoding="utf-8") as file_handle:
                file_handle.write("##fileformat=VCFv4.2\nchr1\tbad\n")

            result = validate_sequence_file(filename)

        self.assertFalse(result["valid"])
        self.assertTrue(any("#CHROM" in error for error in result["errors"]))


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

    def test_pipeline_analyzes_vcf_and_writes_variant_metrics(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "calls.vcf")
            with open(input_file, "w", encoding="utf-8") as file_handle:
                file_handle.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t10\t.\tA\tG\t40\tPASS\t.\n"
                )

            report = run_pipeline(input_file, results_dir=temp_dir)
            self.assertEqual(report["format"], "vcf")
            self.assertEqual(report["metrics"]["variant_records"], 1)
            self.assertEqual(report["metrics"]["snps"], 1)
            self.assertTrue(os.path.isfile(report["report_files"]["json"]))

    def test_pipeline_dispatches_bam_to_alignment_qc(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = os.path.join(temp_dir, "sample.bam")
            with open(input_file, "wb") as file_handle:
                file_handle.write(b"placeholder")

            metrics = {
                "flagstat": {"in total": {"passed": 10, "failed": 0}},
                "stats": {"raw total sequences": 10},
                "reference_file": None,
            }
            with patch("bioseq.samtools.quickcheck"), patch(
                "bioseq.pipeline.analyze_alignment", return_value=metrics
            ) as analyze:
                report = run_pipeline(input_file, results_dir=temp_dir)

        self.assertEqual(report["format"], "bam")
        self.assertEqual(report["metrics"], metrics)
        self.assertIsNone(report["validation"]["record_count"])
        analyze.assert_called_once_with(input_file, reference_file=None)


if __name__ == "__main__":
    unittest.main()
