import argparse
import os

from .fasta import read_fasta, sequence_stats
from .qc import analyze_fastq_quality
from .report import write_report
from .validation import detect_sequence_format, validate_sequence_file
from .variants import analyze_vcf
from .samtools import analyze_alignment


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")


def _resolve_results_dir(results_dir=None):
    if results_dir is None:
        results_dir = DEFAULT_RESULTS_DIR
    elif not os.path.isabs(results_dir):
        results_dir = os.path.join(PROJECT_ROOT, results_dir)
    return os.path.abspath(results_dir)


def run_pipeline(
    input_file,
    results_dir=None,
    minimum_mean_quality=20,
    reference_file=None,
):
    """Validate, analyze, and report FASTA, FASTQ, VCF, BAM, or CRAM input."""
    input_file = os.path.abspath(os.fspath(input_file))
    file_format = detect_sequence_format(input_file)
    if file_format == "vcf":
        vcf_results = analyze_vcf(input_file)
        validation = vcf_results["validation"]
        metrics = vcf_results["metrics"]
    else:
        validation = validate_sequence_file(input_file, file_format=file_format)
    if not validation["valid"]:
        details = "\n".join(f"- {error}" for error in validation["errors"])
        raise ValueError(f"Input validation failed for {input_file}:\n{details}")

    if file_format == "fasta":
        sequences = read_fasta(input_file)
        metrics = sequence_stats(sequences)
        total_bases = sum(len(sequence) for sequence in sequences.values())
        gc_bases = sum(
            sequence.upper().count("G") + sequence.upper().count("C")
            for sequence in sequences.values()
        )
        metrics.update(
            {
                "total_bases": total_bases,
                "gc_content": (gc_bases / total_bases * 100) if total_bases else 0,
            }
        )
    elif file_format == "fastq":
        metrics = analyze_fastq_quality(
            input_file, minimum_mean_quality=minimum_mean_quality
        )
    elif file_format in {"bam", "cram"}:
        metrics = analyze_alignment(input_file, reference_file=reference_file)

    results_dir = _resolve_results_dir(results_dir)
    report_name = f"{os.path.splitext(os.path.basename(input_file))[0]}_report"
    report_data = {
        "input": input_file,
        "format": file_format,
        "reference": os.path.abspath(os.fspath(reference_file))
        if file_format == "cram" and reference_file is not None
        else None,
        "validation": validation,
        "metrics": metrics,
        "report_files": {
            "json": os.path.join(results_dir, f"{report_name}.json"),
            "text": os.path.join(results_dir, f"{report_name}.txt"),
        },
    }
    write_report(report_data, results_dir, report_name)
    return report_data


def main():
    parser = argparse.ArgumentParser(
        description="Validate and analyze FASTA, FASTQ, VCF, BAM, or CRAM files."
    )
    parser.add_argument(
        "input_file", help="FASTA, FASTQ, VCF, BAM, or CRAM file to process"
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help="Directory for JSON and text reports (default: project results/)",
    )
    parser.add_argument(
        "--minimum-mean-quality",
        type=float,
        default=20,
        help="FASTQ mean-read Phred threshold for low-quality reads (default: 20)",
    )
    parser.add_argument(
        "--reference",
        default=None,
        help="Reference FASTA for CRAM decoding; samtools must be installed",
    )
    args = parser.parse_args()

    try:
        report = run_pipeline(
            args.input_file,
            results_dir=args.results_dir,
            minimum_mean_quality=args.minimum_mean_quality,
            reference_file=args.reference,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))

    print(f"Validated and analyzed {report['input']}")
    for report_type, report_path in report["report_files"].items():
        print(f"{report_type.upper()} report: {report_path}")


if __name__ == "__main__":
    main()
