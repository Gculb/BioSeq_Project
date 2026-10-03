import gzip
from io import BytesIO
import os
import tempfile
import unittest
from unittest.mock import patch

from bioseq.download import (
    download_alignment_file,
    download_bam,
    download_fastq,
    download_variant_vcf,
)


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

    def test_download_variant_vcf_saves_and_validates_compressed_vcf(self):
        vcf_data = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\n"
        with tempfile.TemporaryDirectory() as output_dir, patch(
            "bioseq.download.urllib.request.urlopen",
            return_value=BytesIO(gzip.compress(vcf_data)),
        ):
            result = download_variant_vcf(
                "https://example.org/sample.vcf.gz?download=1",
                output_dir=output_dir,
            )

            self.assertEqual(result, os.path.join(output_dir, "sample.vcf.gz"))
            with gzip.open(result, "rb") as file_handle:
                self.assertEqual(file_handle.read(), vcf_data)

    def test_download_bam_saves_and_validates_bam_payload(self):
        bam_data = gzip.compress(b"BAM\x01binary payload")
        with tempfile.TemporaryDirectory() as output_dir, patch(
            "bioseq.download.urllib.request.urlopen",
            return_value=BytesIO(bam_data),
        ):
            result = download_bam(
                "https://example.org/sample.bam", output_dir=output_dir
            )

            self.assertEqual(result, os.path.join(output_dir, "sample.bam"))
            with gzip.open(result, "rb") as file_handle:
                self.assertEqual(file_handle.read(), b"BAM\x01binary payload")

    def test_download_alignment_file_accepts_cram(self):
        with tempfile.TemporaryDirectory() as output_dir, patch(
            "bioseq.download.urllib.request.urlopen",
            return_value=BytesIO(b"CRAM\x03binary payload"),
        ):
            result = download_alignment_file(
                "https://example.org/sample.cram", output_dir=output_dir
            )
        self.assertEqual(result, os.path.join(output_dir, "sample.cram"))

    def test_download_variant_rejects_non_vcf_content(self):
        with tempfile.TemporaryDirectory() as output_dir, patch(
            "bioseq.download.urllib.request.urlopen",
            return_value=BytesIO(b"not a VCF file"),
        ):
            with self.assertRaisesRegex(ValueError, "VCF header"):
                download_variant_vcf(
                    "https://example.org/sample.vcf", output_dir=output_dir
                )
            self.assertEqual(os.listdir(output_dir), [])

    def test_download_rejects_bad_url_and_refuses_to_overwrite(self):
        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(ValueError, "HTTP or HTTPS"):
                download_bam("file:///sample.bam", output_dir=output_dir)

            existing = os.path.join(output_dir, "sample.bam")
            with open(existing, "wb") as file_handle:
                file_handle.write(b"keep this")
            with self.assertRaises(FileExistsError):
                download_bam(
                    "https://example.org/sample.bam", output_dir=output_dir
                )
            with open(existing, "rb") as file_handle:
                self.assertEqual(file_handle.read(), b"keep this")


if __name__ == "__main__":
    unittest.main()
