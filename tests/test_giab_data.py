import os
import tempfile
import unittest

from bioseq.giab_data import _restrict_confident_regions


class GiabRegionPreparationTests(unittest.TestCase):
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
