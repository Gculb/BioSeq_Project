import argparse
import gzip
import hashlib
import json
import os
import tempfile
import urllib.error
import urllib.request

from ._external import require_executable, run_command
from .workflow import _validate_interval


GIAB_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
)
GIAB_READS_BASE = (
    "https://ftp.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/"
    "HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams"
)
TRUTH_VCF_NAME = "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONFIDENT_BED_NAME = "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REFERENCE_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"
)
KNOWN_SITES_URL = (
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/"
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
)
ALIGNMENT_URL = f"{GIAB_READS_BASE}/HG002.GRCh38.2x250.bam"
READ_EXTRACTION_FLANK = 5000


def _download(url, destination):
    request = urllib.request.Request(
        url, headers={"User-Agent": "SeqCheckFlow/1.0 (GIAB benchmark preparation)"}
    )
    temporary_path = f"{destination}.part"
    try:
        with urllib.request.urlopen(request, timeout=120) as response, open(
            temporary_path, "wb"
        ) as output:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
        if os.path.getsize(temporary_path) == 0:
            raise ValueError(f"Downloaded an empty file from {url}")
        os.replace(temporary_path, destination)
    except urllib.error.URLError as exc:
        reason = getattr(exc, "reason", exc)
        raise RuntimeError(f"Failed to download {url}: {reason}") from exc
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
    return destination


def _sha256(filename):
    digest = hashlib.sha256()
    with open(filename, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _restrict_confident_regions(source_bed, output_bed, contig, start, end):
    start_zero_based = start - 1
    end_exclusive = end
    kept = 0
    with open(source_bed, "r", encoding="utf-8") as source, open(
        output_bed, "w", encoding="utf-8"
    ) as output:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3 or fields[0] != contig:
                continue
            region_start = max(int(fields[1]), start_zero_based)
            region_end = min(int(fields[2]), end_exclusive)
            if region_start >= region_end:
                continue
            fields[1] = str(region_start)
            fields[2] = str(region_end)
            output.write("\t".join(fields) + "\n")
            kept += 1
    if kept == 0:
        raise ValueError(
            f"GIAB confident regions do not overlap {contig}:{start}-{end}."
        )
    return kept


def _read_extraction_region(contig, start, end):
    padded_start = max(1, start - READ_EXTRACTION_FLANK)
    padded_end = end + READ_EXTRACTION_FLANK
    return f"{contig}:{padded_start}-{padded_end}"


def _extract_paired_fastq(samtools, threads, regional_bam, read1, read2):
    with tempfile.TemporaryDirectory(prefix="seqcheckflow-giab-") as temporary_directory:
        name_collated_bam = os.path.join(temporary_directory, "regional.name_collated.bam")
        run_command(
            [
                samtools,
                "collate",
                "-@",
                str(threads),
                "-o",
                name_collated_bam,
                regional_bam,
            ]
        )
        run_command(
            [
                samtools,
                "fastq",
                "-@",
                str(threads),
                "-f",
                "1",
                "-F",
                "2304",
                "-n",
                "-1",
                read1,
                "-2",
                read2,
                "-0",
                os.devnull,
                "-s",
                os.devnull,
                name_collated_bam,
            ]
        )


def prepare_giab_region(output_dir, interval="chr20:10000000-11000000", threads=2):
    """Prepare a small paired-read HG002 region from the public GIAB indexed BAM."""
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer.")
    contig, start, end = _validate_interval(interval)
    if contig != "chr20":
        raise ValueError("This sample preparation recipe currently supports chr20 only.")
    output_dir = os.path.abspath(os.fspath(output_dir))
    if os.path.exists(output_dir):
        raise FileExistsError(f"GIAB data directory already exists: {output_dir}")
    os.makedirs(output_dir)

    truth_vcf = _download(
        f"{GIAB_BASE}/{TRUTH_VCF_NAME}",
        os.path.join(output_dir, TRUTH_VCF_NAME),
    )
    _download(
        f"{GIAB_BASE}/{TRUTH_VCF_NAME}.tbi",
        os.path.join(output_dir, f"{TRUTH_VCF_NAME}.tbi"),
    )
    source_bed = _download(
        f"{GIAB_BASE}/{CONFIDENT_BED_NAME}",
        os.path.join(output_dir, CONFIDENT_BED_NAME),
    )
    bed_name = f"HG002_chr20_{start}_{end}.confident.bed"
    confident_bed = os.path.join(output_dir, bed_name)
    confident_interval_count = _restrict_confident_regions(
        source_bed, confident_bed, contig, start, end
    )

    reference_gzip = os.path.join(output_dir, "chr20.fa.gz")
    _download(REFERENCE_URL, reference_gzip)
    reference_fasta = os.path.join(output_dir, "chr20.fa")
    temporary_reference = f"{reference_fasta}.part"
    try:
        with gzip.open(reference_gzip, "rb") as source, open(
            temporary_reference, "wb"
        ) as output:
            while chunk := source.read(1024 * 1024):
                output.write(chunk)
        os.replace(temporary_reference, reference_fasta)
    finally:
        if os.path.exists(temporary_reference):
            os.unlink(temporary_reference)

    samtools = require_executable("samtools")
    bcftools = require_executable("bcftools")
    rtg = require_executable("rtg")
    reference_sdf = os.path.join(output_dir, "chr20.sdf")
    run_command([rtg, "format", "-o", reference_sdf, reference_fasta])

    region = f"{contig}:{start}-{end}"
    regional_truth = os.path.join(output_dir, "HG002.truth.region.vcf.gz")
    run_command(
        [
            bcftools,
            "view",
            "-r",
            region,
            "-Oz",
            "-o",
            regional_truth,
            truth_vcf,
        ]
    )
    run_command([bcftools, "index", "-f", "-t", regional_truth])
    known_sites = os.path.join(output_dir, "Mills.hg38.region.vcf.gz")
    run_command(
        [
            bcftools,
            "view",
            "-r",
            region,
            "-Oz",
            "-o",
            known_sites,
            KNOWN_SITES_URL,
        ]
    )
    run_command([bcftools, "index", "-f", "-t", known_sites])

    regional_bam = os.path.join(output_dir, "HG002.region.bam")
    source_bam = f"{ALIGNMENT_URL}"
    read_extraction_region = _read_extraction_region(contig, start, end)
    run_command(
        [
            samtools,
            "view",
            "-b",
            "-o",
            regional_bam,
            source_bam,
            read_extraction_region,
        ]
    )
    read1 = os.path.join(output_dir, "HG002_R1.fastq")
    read2 = os.path.join(output_dir, "HG002_R2.fastq")
    _extract_paired_fastq(samtools, threads, regional_bam, read1, read2)
    if not os.path.getsize(read1) or not os.path.getsize(read2):
        raise RuntimeError("GIAB region extraction produced an empty FASTQ mate.")

    return {
        "sample": "HG002",
        "region": region,
        "read_extraction_region": read_extraction_region,
        "read1": read1,
        "read2": read2,
        "reference_fasta": reference_fasta,
        "reference_sdf": reference_sdf,
        "known_sites_vcf": known_sites,
        "truth_vcf": regional_truth,
        "confident_regions_bed": confident_bed,
        "source_alignment_url": source_bam,
        "source_truth_vcf_url": f"{GIAB_BASE}/{TRUTH_VCF_NAME}",
        "source_truth_bed_url": f"{GIAB_BASE}/{CONFIDENT_BED_NAME}",
        "source_reference_url": REFERENCE_URL,
        "source_known_sites_url": KNOWN_SITES_URL,
        "source_alignment_access": (
            "The indexed 122-GB public BAM is queried by genomic interval with "
            "5-kb flanks; the whole BAM is not intentionally downloaded. "
            "Mate fetching across the full BAM is avoided; the regional BAM is "
            "name-collated before paired FASTQ extraction."
        ),
        "confident_interval_count": confident_interval_count,
        "sha256": {
            "read1": _sha256(read1),
            "read2": _sha256(read2),
            "reference_fasta": _sha256(reference_fasta),
            "known_sites_vcf": _sha256(known_sites),
            "truth_vcf": _sha256(regional_truth),
            "confident_regions_bed": _sha256(confident_bed),
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description="Prepare a region-sized public GIAB HG002 benchmark dataset."
    )
    parser.add_argument(
        "--interval", default="chr20:10000000-11000000", help="GRCh38 interval"
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--output-dir", default="data/benchmark/giab_hg002_chr20"
    )
    args = parser.parse_args()
    try:
        report = prepare_giab_region(
            args.output_dir, interval=args.interval, threads=args.threads
        )
    except (FileExistsError, OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    report_path = os.path.join(os.path.abspath(args.output_dir), "dataset.json")
    with open(report_path, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=2, sort_keys=True)
        file_handle.write("\n")
    print(f"Prepared {report['sample']} {report['region']}: {report_path}")


if __name__ == "__main__":
    main()
