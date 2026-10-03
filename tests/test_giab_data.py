import os
import tempfile
import unittest
from unittest.mock import patch

from bioseq.giab_data import (
    _extract_paired_fastq,
    _read_extraction_region,
    _restrict_confident_regions,
)


class GiabRegionPreparationTests(unittest.TestCase):
    def test_collates_regional_alignments_before_writing_paired_fastqs(self):
        with tempfile.TemporaryDirectory() as directory:
            read1 = os.path.join(directory, "read1.fastq")
            read2 = os.path.join(directory, "read2.fastq")

            def mock_run(command):
                if command[1] == "collate":
                    with open(command[command.index("-o") + 1], "wb") as bam:
                        bam.write(b"collated")
                else:
                    for path_flag in ("-1", "-2"):
                        with open(
                            command[command.index(path_flag) + 1],
                            "w",
                            encoding="utf-8",
                        ) as fastq:
                            fastq.write("@read\nACGT\n+\nIIII\n")
                return ""

            with patch("bioseq.giab_data.run_command", side_effect=mock_run) as run:
                _extract_paired_fastq("samtools", 2, "regional.bam", read1, read2)

            self.assertEqual(run.call_count, 2)
            collate_command = run.call_args_list[0].args[0]
            fastq_command = run.call_args_list[1].args[0]
            self.assertEqual(collate_command[1], "collate")
            self.assertEqual(fastq_command[1], "fastq")
            self.assertEqual(fastq_command[-1], collate_command[collate_command.index("-o") + 1])
            self.assertTrue(os.path.isfile(read1))
            self.assertTrue(os.path.isfile(read2))

    def test_read_extraction_region_includes_flanks_without_negative_start(self):
        self.assertEqual(
            _read_extraction_region("chr20", 10000000, 11000000),
            "chr20:9995000-11005000",
        )
        self.assertEqual(
            _read_extraction_region("chr20", 100, 200),
            "chr20:1-5200",
        )

    def test_restricts_and_clips_bed_intervals_to_requested_region(self):
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.bed")
            target = os.path.join(directory, "target.bed")
            with open(source, "w", encoding="utf-8") as file_handle:
                file_handle.write(
                    "chr19\t0\t500\n"
                    "chr20\t9999990\t10000010\n"
                    "chr20\t10000020\t10000030\tregionA\n"
                    "chr20\t11000000\t11000010\n"
                )

            count = _restrict_confident_regions(
                source, target, "chr20", 10000000, 11000000
            )

            self.assertEqual(count, 2)
            with open(target, encoding="utf-8") as file_handle:
                lines = file_handle.read().splitlines()
        self.assertEqual(lines[0], "chr20\t9999999\t10000010")
        self.assertEqual(lines[1], "chr20\t10000020\t10000030\tregionA")

    def test_rejects_region_without_truth_confidence(self):
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "source.bed")
            target = os.path.join(directory, "target.bed")
            with open(source, "w", encoding="utf-8") as file_handle:
                file_handle.write("chr1\t0\t100\n")

            with self.assertRaisesRegex(ValueError, "do not overlap"):
                _restrict_confident_regions(source, target, "chr20", 1, 100)


if __name__ == "__main__":
    unittest.main()
