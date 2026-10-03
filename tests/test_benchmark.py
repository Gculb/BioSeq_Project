import json
import os
import tempfile
import unittest
from unittest.mock import patch

from bioseq.benchmark import benchmark_variant_calls


SUMMARY = (
    "Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure\n"
    "----------------------------------------------------------------------------------------------------\n"
    "None       8                  8              2          1          0.8000     0.8889       0.8421\n"
)


def write_vcf(filename):
    with open(filename, "w", encoding="utf-8") as file_handle:
        file_handle.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "chr1\t10\t.\tA\tG\t40\tPASS\t.\tGT\t0/1\n"
        )


class VariantBenchmarkTests(unittest.TestCase):
    def test_compares_candidate_and_baseline_and_writes_report(self):
        with tempfile.TemporaryDirectory() as directory:
            truth = os.path.join(directory, "truth.vcf")
            baseline = os.path.join(directory, "baseline.vcf")
            candidate = os.path.join(directory, "candidate.vcf")
            reference = os.path.join(directory, "reference.sdf")
            output_dir = os.path.join(directory, "benchmark")
            os.makedirs(reference)
            for filename in (truth, baseline, candidate):
                write_vcf(filename)

            def run_rtg(command):
                output = command[command.index("--output") + 1]
                os.makedirs(output)
                summary = SUMMARY
                if output.endswith("candidate"):
                    summary = summary.replace(
                        "0.8000     0.8889       0.8421",
                        "0.9000     0.8889       0.8947",
                    )
                with open(
                    os.path.join(output, "summary.txt"), "w", encoding="utf-8"
                ) as file_handle:
                    file_handle.write(summary)

            with patch(
                "bioseq.benchmark.require_executable", return_value="rtg"
            ), patch("bioseq.benchmark.run_command", side_effect=run_rtg) as run_command:
                report = benchmark_variant_calls(
                    truth, baseline, candidate, reference, output_dir
                )

            self.assertEqual(run_command.call_count, 2)
            self.assertEqual(report["baseline"]["metrics"]["false_positives"], 2)
            self.assertAlmostEqual(
                report["candidate"]["metrics"]["precision"], 0.9
            )
            self.assertAlmostEqual(
                report["delta_candidate_minus_baseline"]["precision"], 0.1
            )
            with open(report["report_file"], encoding="utf-8") as file_handle:
                saved_report = json.load(file_handle)
            self.assertEqual(saved_report["candidate"]["calls_vcf"], candidate)

    def test_passes_confident_regions_to_both_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            truth = os.path.join(directory, "truth.vcf")
            baseline = os.path.join(directory, "baseline.vcf")
            candidate = os.path.join(directory, "candidate.vcf")
            regions = os.path.join(directory, "confident.bed")
            reference = os.path.join(directory, "reference.sdf")
            output_dir = os.path.join(directory, "benchmark")
            os.makedirs(reference)
            for filename in (truth, baseline, candidate):
                write_vcf(filename)
            with open(regions, "w", encoding="utf-8") as file_handle:
                file_handle.write("chr1\t0\t20\n")

            def run_rtg(command):
                output = command[command.index("--output") + 1]
                os.makedirs(output)
                with open(
                    os.path.join(output, "summary.txt"), "w", encoding="utf-8"
                ) as file_handle:
                    file_handle.write(SUMMARY)

            with patch(
                "bioseq.benchmark.require_executable", return_value="rtg"
            ), patch("bioseq.benchmark.run_command", side_effect=run_rtg) as run_command:
                benchmark_variant_calls(
                    truth,
                    baseline,
                    candidate,
                    reference,
                    output_dir,
                    regions_bed=regions,
                )

            for call in run_command.call_args_list:
                command = call.args[0]
                self.assertIn("--evaluation-regions", command)
                self.assertIn(regions, command)

    def test_rejects_existing_benchmark_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            truth = os.path.join(directory, "truth.vcf")
            baseline = os.path.join(directory, "baseline.vcf")
            candidate = os.path.join(directory, "candidate.vcf")
            reference = os.path.join(directory, "reference.sdf")
            output_dir = os.path.join(directory, "benchmark")
            os.makedirs(reference)
            os.makedirs(os.path.join(output_dir, "baseline"))
            for filename in (truth, baseline, candidate):
                write_vcf(filename)

            with self.assertRaisesRegex(FileExistsError, "already exists"):
                benchmark_variant_calls(
                    truth, baseline, candidate, reference, output_dir
                )

    def test_rejects_invalid_summary_without_success_shaped_report(self):
        with tempfile.TemporaryDirectory() as directory:
            truth = os.path.join(directory, "truth.vcf")
            baseline = os.path.join(directory, "baseline.vcf")
            candidate = os.path.join(directory, "candidate.vcf")
            reference = os.path.join(directory, "reference.sdf")
            output_dir = os.path.join(directory, "benchmark")
            os.makedirs(reference)
            for filename in (truth, baseline, candidate):
                write_vcf(filename)

            def run_rtg(command):
                output = command[command.index("--output") + 1]
                os.makedirs(output)
                with open(
                    os.path.join(output, "summary.txt"), "w", encoding="utf-8"
                ) as file_handle:
                    file_handle.write(
                        "Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure\n"
                        "----------------------------------------------------------------------------------------------------\n"
                    )

            with patch(
                "bioseq.benchmark.require_executable", return_value="rtg"
            ), patch("bioseq.benchmark.run_command", side_effect=run_rtg):
                with self.assertRaisesRegex(RuntimeError, "no unthresholded metric row"):
                    benchmark_variant_calls(
                        truth, baseline, candidate, reference, output_dir
                    )

            self.assertFalse(os.path.exists(os.path.join(output_dir, "benchmark.json")))

    def test_requires_single_sample_vcfs(self):
        with tempfile.TemporaryDirectory() as directory:
            truth = os.path.join(directory, "truth.vcf")
            baseline = os.path.join(directory, "baseline.vcf")
            candidate = os.path.join(directory, "candidate.vcf")
            reference = os.path.join(directory, "reference.sdf")
            output_dir = os.path.join(directory, "benchmark")
            os.makedirs(reference)
            for filename in (truth, baseline, candidate):
                write_vcf(filename)
            with open(truth, "w", encoding="utf-8") as file_handle:
                file_handle.write(
                    "##fileformat=VCFv4.2\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n"
                    "chr1\t10\t.\tA\tG\t40\tPASS\t.\tGT\t0/1\t0/1\n"
                )

            with self.assertRaisesRegex(ValueError, "exactly one sample"):
                benchmark_variant_calls(
                    truth, baseline, candidate, reference, output_dir
                )


if __name__ == "__main__":
    unittest.main()
