import json
import os
import tempfile
import unittest
from unittest.mock import patch

from bioseq.full_benchmark import compare_germline_workflows


class FullBenchmarkTests(unittest.TestCase):
    def test_runs_baseline_and_bioseq_then_writes_combined_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            read1 = os.path.join(directory, "read1.fastq")
            read2 = os.path.join(directory, "read2.fastq")
            reference = os.path.join(directory, "reference.fa")
            known_sites = os.path.join(directory, "known_sites.vcf.gz")
            truth = os.path.join(directory, "truth.vcf")
            bed = os.path.join(directory, "regions.bed")
            reference_sdf = os.path.join(directory, "reference.sdf")
            output = os.path.join(directory, "comparison")
            for path in (read1, read2, reference, known_sites, truth, bed):
                with open(path, "w", encoding="utf-8") as file_handle:
                    file_handle.write("fixture\n")
            os.makedirs(reference_sdf)
            with open(os.path.join(reference_sdf, "reference.sdf"), "w") as file_handle:
                file_handle.write("fixture\n")
            workflow_results = [
                {
                    "implementation": "baseline",
                    "report_file": "baseline/workflow.json",
                    "variant_vcf": "baseline/calls.vcf.gz",
                    "total_wall_seconds": 12.0,
                    "peak_stage_rss_bytes": 1000,
                    "total_output_bytes": 2000,
                    "tool_versions": {"bwa": "0.7.17", "samtools": "1.20"},
                    "stages": [
                        {
                            "name": "bwa",
                            "wall_seconds": 4.0,
                            "peak_process_tree_rss_bytes": 800,
                            "process_tree_cpu_seconds": 3.5,
                            "output_bytes": 1000,
                        }
                    ],
                },
                {
                    "implementation": "bioseq",
                    "report_file": "bioseq/workflow.json",
                    "variant_vcf": "bioseq/calls.vcf.gz",
                    "total_wall_seconds": 13.0,
                    "peak_stage_rss_bytes": 1100,
                    "total_output_bytes": 2000,
                    "tool_versions": {"bwa": "0.7.17", "samtools": "1.20"},
                    "stages": [
                        {
                            "name": "bwa",
                            "wall_seconds": 4.5,
                            "peak_process_tree_rss_bytes": 900,
                            "process_tree_cpu_seconds": 3.6,
                            "output_bytes": 1000,
                        }
                    ],
                },
            ]
            evaluation = {
                "baseline": {"metrics": {"f1_score": 0.9}},
                "candidate": {"metrics": {"f1_score": 0.9}},
            }
            with patch(
                "bioseq.full_benchmark.run_germline_workflow",
                side_effect=workflow_results,
            ) as run_workflow, patch(
                "bioseq.full_benchmark.benchmark_variant_calls",
                return_value=evaluation,
            ):
                result = compare_germline_workflows(
                    read1,
                    read2,
                    reference,
                    known_sites,
                    reference_sdf,
                    "chr20:10000000-11000000",
                    truth,
                    bed,
                    output,
                )

            self.assertEqual(run_workflow.call_count, 2)
            self.assertEqual(
                result["candidate_minus_baseline"]["total_wall_seconds"], 1.0
            )
            self.assertEqual(
                result["candidate_minus_baseline"]["stages"]["bwa"]["wall_seconds"],
                0.5,
            )
            self.assertEqual(result["truth_evaluation"], evaluation)
            with open(result["report_file"], encoding="utf-8") as file_handle:
                self.assertEqual(json.load(file_handle)["region"], "chr20:10000000-11000000")


if __name__ == "__main__":
    unittest.main()
