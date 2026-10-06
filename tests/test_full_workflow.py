import json
import os
import tempfile
import unittest
from unittest.mock import patch

from seqcheckflow.workflow import run_germline_workflow


def make_fastq(filename, sequence="ACGT"):
    with open(filename, "w", encoding="utf-8") as file_handle:
        file_handle.write(f"@read\n{sequence}\n+\nIIII\n")


def make_file(filename, content=b"fixture"):
    with open(filename, "wb") as file_handle:
        file_handle.write(content)


class GermlineWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.read1 = os.path.join(self.temp_dir.name, "read1.fastq")
        self.read2 = os.path.join(self.temp_dir.name, "read2.fastq")
        self.reference = os.path.join(self.temp_dir.name, "reference.fa")
        self.known_sites = os.path.join(self.temp_dir.name, "known_sites.vcf.gz")
        make_fastq(self.read1)
        make_fastq(self.read2, "TGCA")
        make_file(self.reference, b">chr20\nACGT\n")
        make_file(self.known_sites, b"known sites")

    def _mock_common_tools(self):
        def run(command, stdout_path=None):
            if stdout_path:
                make_file(stdout_path)
            elif command[0] == "bwa" and command[1] == "index":
                make_file(f"{command[-1]}.bwt")
            elif "faidx" in command:
                make_file(f"{command[-1]}.fai")
            elif command[0] == "samtools" and command[1] == "fixmate":
                make_file(command[-1])
            elif command[0] == "samtools" and "-o" in command:
                make_file(command[command.index("-o") + 1])
            elif "-O" in command:
                make_file(command[command.index("-O") + 1])
            elif "-o" in command:
                make_file(command[command.index("-o") + 1])
            else:
                make_file(command[-1])
            return ""

        return patch("seqcheckflow.workflow.run_command", side_effect=run)

    def test_seqcheckflow_workflow_runs_wrappers_and_records_stage_metrics(self):
        output_dir = os.path.join(self.temp_dir.name, "seqcheckflow")
        metrics = {"variant_records": 1, "snps": 1}
        vcf_result = {
            "validation": {
                "valid": True,
                "format": "vcf",
                "record_count": 1,
                "errors": [],
            },
            "metrics": metrics,
        }
        with patch(
            "seqcheckflow.workflow.require_executable", side_effect=lambda name: name
        ), patch("seqcheckflow.workflow.index_reference") as index_ref, patch(
            "seqcheckflow.workflow.align_reads",
            side_effect=lambda reads, ref, output, **kwargs: (
                make_file(output),
                output,
            )[1],
        ) as align, patch(
            "seqcheckflow.workflow.convert_sam_to_bam",
            side_effect=lambda source, output: (make_file(output), output)[1],
        ), patch(
            "seqcheckflow.workflow.sort_bam",
            side_effect=lambda source, output, by_name=False: (
                make_file(output),
                output,
            )[1],
        ), patch(
            "seqcheckflow.workflow.fixmate",
            side_effect=lambda source, output: (make_file(output), output)[1],
        ), patch(
            "seqcheckflow.workflow.mark_duplicates",
            side_effect=lambda source, output: (make_file(output), output)[1],
        ), patch(
            "seqcheckflow.workflow.index_bam",
            side_effect=lambda bam: (make_file(f"{bam}.bai"), f"{bam}.bai")[1],
        ), patch(
            "seqcheckflow.workflow.analyze_alignment",
            return_value={"flagstat": {"mapped": 1}},
        ), patch(
            "seqcheckflow.workflow._run_gatk_haplotype_caller",
            side_effect=lambda ref, bam, interval, threads, output: (
                make_file(output),
                output,
            )[1],
        ), patch(
            "seqcheckflow.workflow.analyze_vcf", return_value=vcf_result
        ), patch(
            "seqcheckflow.workflow._tool_versions",
            return_value={
                "bwa": "test",
                "samtools": "test",
                "gatk": "test",
                "rtg": "test",
            },
        ), self._mock_common_tools():
            result = run_germline_workflow(
                self.read1,
                self.read2,
                self.reference,
                self.known_sites,
                "chr20:10000000-11000000",
                output_dir,
            )

        self.assertEqual(result["implementation"], "seqcheckflow")
        self.assertEqual(result["variant_metrics"], metrics)
        self.assertEqual(
            result["inputs"]["reference"]["size_bytes"], os.path.getsize(self.reference)
        )
        self.assertIn(
            "bwa_mem_alignment", [stage["name"] for stage in result["stages"]]
        )
        self.assertTrue(os.path.isfile(result["report_file"]))
        with open(result["report_file"], encoding="utf-8") as file_handle:
            self.assertEqual(json.load(file_handle)["sample_name"], "HG002")
        self.assertTrue(align.call_args.kwargs["read_group"].startswith("@RG\\t"))
        index_ref.assert_called_once()

    def test_direct_tool_baseline_runs_bqsr_and_gatk_without_seqcheckflow_wrappers(self):
        output_dir = os.path.join(self.temp_dir.name, "baseline")
        vcf_result = {
            "validation": {
                "valid": True,
                "format": "vcf",
                "record_count": 1,
                "errors": [],
            },
            "metrics": {"variant_records": 1},
        }
        with patch(
            "seqcheckflow.workflow.require_executable", side_effect=lambda name: name
        ), patch(
            "seqcheckflow.workflow._tool_versions",
            return_value={
                "bwa": "test",
                "samtools": "test",
                "gatk": "test",
                "rtg": "test",
            },
        ), patch(
            "seqcheckflow.workflow.analyze_alignment", return_value={"coverage": []}
        ), patch(
            "seqcheckflow.workflow._run_gatk_haplotype_caller",
            side_effect=lambda ref, bam, interval, threads, output: (
                make_file(output),
                output,
            )[1],
        ), patch(
            "seqcheckflow.workflow.analyze_vcf", return_value=vcf_result
        ), self._mock_common_tools() as run_command:
            result = run_germline_workflow(
                self.read1,
                self.read2,
                self.reference,
                self.known_sites,
                "chr20:10000000-11000000",
                output_dir,
                implementation="baseline",
            )

        commands = [call.args[0] for call in run_command.call_args_list]
        self.assertTrue(
            any(command[:2] == ["bwa", "mem"] for command in commands)
        )
        self.assertTrue(
            any(command[:2] == ["samtools", "markdup"] for command in commands)
        )
        self.assertTrue(
            any(command[1] == "BaseRecalibrator" for command in commands)
        )
        self.assertTrue(any(command[1] == "ApplyBQSR" for command in commands))
        self.assertTrue(result["variant_vcf"].endswith("calls.vcf.gz"))

    def test_rejects_invalid_interval_and_existing_output(self):
        with self.assertRaisesRegex(ValueError, "interval"):
            run_germline_workflow(
                self.read1,
                self.read2,
                self.reference,
                self.known_sites,
                "chr20:0-10",
                os.path.join(self.temp_dir.name, "invalid"),
            )

        existing = os.path.join(self.temp_dir.name, "existing")
        os.makedirs(existing)
        with self.assertRaisesRegex(FileExistsError, "already exists"):
            run_germline_workflow(
                self.read1,
                self.read2,
                self.reference,
                self.known_sites,
                "chr20:10000000-11000000",
                existing,
            )


if __name__ == "__main__":
    unittest.main()
