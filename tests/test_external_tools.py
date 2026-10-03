from io import StringIO
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from bioseq.alignment import align_reads, convert_sam_to_bam
from bioseq.blast import BLAST_search
from bioseq.samtools import (
    analyze_alignment,
    depth,
    flagstat,
    index_bam,
    sort_bam,
    stats,
)


def successful_tool(command):
    output_path = command[-2]
    with open(output_path, "wb") as file_handle:
        file_handle.write(b"tool output")
    return SimpleNamespace(stdout="", stderr="")


class BlastTests(unittest.TestCase):
    def test_blast_returns_one_summary_per_hit(self):
        best_hsp = SimpleNamespace(
            expect=1e-20,
            align_length=40,
            identities=36,
            query_start=11,
            query_end=50,
            bits=85.5,
        )
        record = SimpleNamespace(
            query_length=100,
            alignments=[
                SimpleNamespace(
                    hit_id="gi|123",
                    accession="ABC123",
                    hit_def="example hit",
                    hsps=[best_hsp],
                )
            ],
        )
        response = StringIO("<BlastOutput/>")
        with patch("bioseq.blast.NCBIWWW.qblast", return_value=response) as qblast, patch(
            "bioseq.blast.NCBIXML.read", return_value=record
        ):
            hits = BLAST_search(
                "ACGT",
                email="researcher@example.org",
                hitlist_size=5,
                expect=0.01,
            )

        self.assertEqual(hits[0]["hit_id"], "gi|123")
        self.assertEqual(hits[0]["accession"], "ABC123")
        self.assertEqual(hits[0]["percent_identity"], 90)
        self.assertEqual(hits[0]["query_coverage"], 40)
        self.assertEqual(hits[0]["evalue"], 1e-20)
        self.assertTrue(response.closed)
        self.assertEqual(qblast.call_args.args, ("blastn", "nt", "ACGT"))
        self.assertEqual(qblast.call_args.kwargs["email"], "researcher@example.org")

    def test_blast_validates_query_and_program(self):
        with self.assertRaisesRegex(ValueError, "non-empty"):
            BLAST_search("  ")
        with self.assertRaisesRegex(ValueError, "Unsupported BLAST program"):
            BLAST_search("ACGT", program="blastz")


class AlignmentTests(unittest.TestCase):
    def test_align_reads_builds_bwa_mem_command(self):
        with tempfile.TemporaryDirectory() as directory:
            reads = os.path.join(directory, "reads.fastq")
            reference = os.path.join(directory, "reference.fasta")
            output = os.path.join(directory, "aligned.sam")
            for path in (reads, reference, f"{reference}.bwt"):
                with open(path, "w", encoding="utf-8") as file_handle:
                    file_handle.write("fixture\n")

            with patch("bioseq.alignment.require_executable", return_value="bwa"), patch(
                "bioseq.alignment.run_command", return_value=output
            ) as run_command:
                result = align_reads(reads, reference, output, threads=3)

        self.assertEqual(result, output)
        self.assertEqual(
            run_command.call_args.args[0],
            ["bwa", "mem", "-t", "3", reference, reads],
        )
        self.assertEqual(run_command.call_args.kwargs["stdout_path"], output)

    def test_bwa_requires_an_index(self):
        with tempfile.TemporaryDirectory() as directory:
            reads = os.path.join(directory, "reads.fastq")
            reference = os.path.join(directory, "reference.fasta")
            for path in (reads, reference):
                with open(path, "w", encoding="utf-8") as file_handle:
                    file_handle.write("fixture\n")

            with patch("bioseq.alignment.require_executable", return_value="bwa"):
                with self.assertRaisesRegex(FileNotFoundError, "bwa index"):
                    align_reads(reads, reference, os.path.join(directory, "out.sam"))

    def test_align_reads_supports_paired_end_files(self):
        with tempfile.TemporaryDirectory() as directory:
            reads = [os.path.join(directory, f"reads_{mate}.fastq") for mate in (1, 2)]
            reference = os.path.join(directory, "reference.fasta")
            output = os.path.join(directory, "aligned.sam")
            for path in [*reads, reference]:
                with open(path, "w", encoding="utf-8") as file_handle:
                    file_handle.write("fixture\n")

            with patch(
                "bioseq.alignment.require_executable", return_value="minimap2"
            ), patch(
                "bioseq.alignment.run_command", return_value=output
            ) as run_command:
                result = align_reads(
                    reads, reference, output, aligner="minimap2", threads=2
                )

        self.assertEqual(result, output)
        self.assertEqual(
            run_command.call_args.args[0],
            ["minimap2", "-a", "-t", "2", reference, *reads],
        )

    def test_sam_conversion_writes_and_returns_bam_path(self):
        with tempfile.TemporaryDirectory() as directory:
            sam = os.path.join(directory, "alignments.sam")
            with open(sam, "w", encoding="utf-8") as file_handle:
                file_handle.write("@HD\n")

            with patch("bioseq.alignment.require_executable", return_value="samtools"), patch(
                "bioseq.alignment.run_command", side_effect=successful_tool
            ) as run_command:
                result = convert_sam_to_bam(sam)

            self.assertEqual(result, os.path.join(directory, "alignments.bam"))
            self.assertTrue(os.path.isfile(result))
            self.assertEqual(run_command.call_args.args[0][:3], ["samtools", "view", "-b"])


class SamtoolsTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.bam = os.path.join(self.temp_dir.name, "alignments.bam")
        with open(self.bam, "wb") as file_handle:
            file_handle.write(b"bam fixture")
        self.tool_patch = patch("bioseq.samtools.require_executable", return_value="samtools")
        self.tool_patch.start()
        self.addCleanup(self.tool_patch.stop)
        self.addCleanup(self.temp_dir.cleanup)

    def test_sort_bam_creates_sorted_output(self):
        with patch("bioseq.samtools.run_command", side_effect=successful_tool):
            sorted_bam = sort_bam(self.bam)

        self.assertEqual(sorted_bam, os.path.join(self.temp_dir.name, "alignments.sorted.bam"))
        self.assertTrue(os.path.isfile(sorted_bam))

    def test_index_bam_creates_index(self):
        with patch("bioseq.samtools.run_command", side_effect=successful_tool):
            index_file = index_bam(self.bam)

        self.assertEqual(index_file, f"{self.bam}.bai")
        self.assertTrue(os.path.isfile(index_file))

    def test_flagstat_parses_counts_and_percentages(self):
        output = (
            "12 + 1 in total (QC-passed reads + QC-failed reads)\n"
            "10 + 1 mapped (83.33% : 100.00%)\n"
            "2 + 0 unmapped\n"
        )
        with patch("bioseq.samtools.run_command", return_value=output):
            result = flagstat(self.bam)

        self.assertEqual(result["in total"]["passed"], 12)
        self.assertEqual(result["in total"]["failed"], 1)
        self.assertEqual(result["mapped"]["percentages"], [83.33, 100.0])
        self.assertEqual(result["unmapped"]["passed"], 2)

    def test_stats_parses_summary_values(self):
        output = (
            "SN\traw total sequences:\t12\n"
            "SN\treads mapped:\t10\n"
            "SN\taverage length:\t150.5\n"
            "LF\tignored format line\n"
        )
        with patch("bioseq.samtools.run_command", return_value=output):
            result = stats(self.bam)

        self.assertEqual(result, {
            "raw total sequences": 12,
            "reads mapped": 10,
            "average length": 150.5,
        })

    def test_depth_can_return_text_or_write_a_file(self):
        depth_text = "chr1\t1\t4\nchr1\t2\t8\n"
        with patch("bioseq.samtools.run_command", return_value=depth_text) as run_command:
            result = depth(self.bam)
        self.assertEqual(result, depth_text)
        self.assertEqual(run_command.call_args.args[0][:2], ["samtools", "depth"])

        output_path = os.path.join(self.temp_dir.name, "depth.tsv")
        def write_depth(command, stdout_path):
            self.assertEqual(command[:2], ["samtools", "depth"])
            with open(stdout_path, "w", encoding="utf-8") as file_handle:
                file_handle.write(depth_text)
            return stdout_path

        with patch("bioseq.samtools.run_command", side_effect=write_depth):
            result_path = depth(self.bam, output_file=output_path)
        self.assertEqual(result_path, output_path)
        with open(output_path, encoding="utf-8") as file_handle:
            self.assertEqual(file_handle.read(), depth_text)

    def test_analyze_cram_uses_supplied_reference_for_temporary_bam(self):
        cram_file = os.path.join(self.temp_dir.name, "alignments.cram")
        reference_file = os.path.join(self.temp_dir.name, "reference.fasta")
        with open(cram_file, "wb") as file_handle:
            file_handle.write(b"CRAM")
        with open(reference_file, "w", encoding="utf-8") as file_handle:
            file_handle.write(">chr1\nACGT\n")

        def convert_cram(command):
            self.assertEqual(command[1:5], ["view", "-T", reference_file, "-b"])
            output_path = command[6]
            with open(output_path, "wb") as file_handle:
                file_handle.write(b"decoded bam")
            return ""

        with patch("bioseq.samtools.quickcheck"), patch(
            "bioseq.samtools.run_command", side_effect=convert_cram
        ), patch(
            "bioseq.samtools.flagstat", return_value={"mapped": {"passed": 8}}
        ) as mock_flagstat, patch(
            "bioseq.samtools.stats", return_value={"raw total sequences": 10}
        ) as mock_stats:
            result = analyze_alignment(cram_file, reference_file=reference_file)

        self.assertEqual(result["flagstat"], {"mapped": {"passed": 8}})
        self.assertEqual(result["stats"], {"raw total sequences": 10})
        self.assertEqual(result["reference_file"], reference_file)
        analyzed_file = mock_flagstat.call_args.args[0]
        self.assertEqual(analyzed_file, mock_stats.call_args.args[0])
        self.assertTrue(analyzed_file.endswith("decoded.bam"))


if __name__ == "__main__":
    unittest.main()
