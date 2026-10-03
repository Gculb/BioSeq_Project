import os
import tempfile
import unittest
from unittest.mock import patch

from bioseq.download import download_fastq


class DownloadFastqTests(unittest.TestCase):
    def test_download_fastq_returns_sra_output_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            with open(os.path.join(output_dir, "SRR1234567.fastq"), "w") as handle:
                handle.write("existing file with a similar prefix\n")

            def run_fasterq(command, **kwargs):
                self.assertIn("--split-files", command)
                self.assertTrue(kwargs["check"])
                with open(os.path.join(output_dir, "SRR123456_1.fastq"), "w") as handle:
                    handle.write("@read/1\nACGT\n+\nIIII\n")
                with open(os.path.join(output_dir, "SRR123456_2.fastq"), "w") as handle:
                    handle.write("@read/2\nTGCA\n+\nIIII\n")

            with patch("bioseq.download.shutil.which", return_value="fasterq-dump"), patch(
                "bioseq.download.subprocess.run", side_effect=run_fasterq
            ) as run:
                paths = download_fastq("SRR123456", output_dir=output_dir, threads=2)

        self.assertEqual(
            paths,
            [
                os.path.join(output_dir, "SRR123456_1.fastq"),
                os.path.join(output_dir, "SRR123456_2.fastq"),
            ],
        )
        self.assertIn("--split-files", run.call_args.args[0])
        self.assertIn("2", run.call_args.args[0])
        self.assertTrue(run.call_args.kwargs["check"])

    def test_download_fastq_rejects_non_run_accession(self):
        with self.assertRaisesRegex(ValueError, "SRA run accession"):
            download_fastq("NC_000913")

    def test_download_fastq_reports_missing_sra_toolkit(self):
        with patch("bioseq.download.shutil.which", return_value=None):
            with self.assertRaisesRegex(FileNotFoundError, "SRA Toolkit"):
                download_fastq("SRR123456")


if __name__ == "__main__":
    unittest.main()
