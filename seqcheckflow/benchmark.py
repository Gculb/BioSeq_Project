import argparse
import json
import math
import os

from ._external import require_executable, require_file, run_command
from .variants import analyze_vcf


SUMMARY_METRICS = {
    "true-pos-baseline": ("true_positives_baseline", int),
    "true-pos-call": ("true_positives_candidate", int),
    "false-pos": ("false_positives", int),
    "false-neg": ("false_negatives", int),
    "precision": ("precision", float),
    "sensitivity": ("sensitivity", float),
    "f-measure": ("f1_score", float),
}


def _validate_vcf(filename, description):
    filename = require_file(filename, description)
    result = analyze_vcf(filename)
    if not result["validation"]["valid"]:
        details = "\n".join(f"- {error}" for error in result["validation"]["errors"])
        raise ValueError(f"{description} VCF is invalid:\n{details}")
    if result["metrics"]["sample_count"] != 1:
        raise ValueError(
            f"{description} VCF must contain exactly one sample for this benchmark."
        )
    return filename


def _parse_vcfeval_summary(filename):
    headers = None
    values = None
    with open(filename, "r", encoding="utf-8") as file_handle:
        for line in file_handle:
            fields = line.split()
            normalized_fields = [field.casefold() for field in fields]
            if all(key in normalized_fields for key in SUMMARY_METRICS):
                headers = normalized_fields
                continue
            if headers is None or not fields or fields[0].casefold() != "none":
                continue
            if len(fields) != len(headers):
                raise RuntimeError(f"Invalid RTG summary row in {filename}: {line.rstrip()}")
            values = dict(zip(headers, fields))
            break

    if values is None:
        raise RuntimeError(
            f"RTG summary {filename} has no unthresholded metric row."
        )
    metrics = {}
    for key, (metric_name, value_type) in SUMMARY_METRICS.items():
        try:
            value = value_type(values[key])
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid {key} value in RTG summary {filename}: {values[key]!r}"
            ) from exc
        if isinstance(value, float) and not math.isfinite(value):
            raise RuntimeError(f"Non-finite {key} value in RTG summary {filename}.")
        metrics[metric_name] = value
    return metrics


def _run_vcfeval(
    executable,
    truth_vcf,
    calls_vcf,
    reference_sdf,
    output_dir,
    regions_bed=None,
):
    command = [
        executable,
        "vcfeval",
        "--baseline",
        truth_vcf,
        "--calls",
        calls_vcf,
        "--template",
        reference_sdf,
        "--output",
        output_dir,
        "--no-roc",
    ]
    if regions_bed is not None:
        command.extend(["--evaluation-regions", regions_bed])

    run_command(command)
    summary_file = os.path.join(output_dir, "summary.txt")
    if not os.path.isfile(summary_file):
        raise RuntimeError(f"RTG vcfeval did not create its summary: {summary_file}")
    return _parse_vcfeval_summary(summary_file), summary_file


def benchmark_variant_calls(
    truth_vcf,
    baseline_vcf,
    candidate_vcf,
    reference_sdf,
    output_dir,
    regions_bed=None,
):
    """Compare baseline and candidate VCFs against truth using RTG vcfeval."""
    truth_vcf = _validate_vcf(truth_vcf, "Truth")
    baseline_vcf = _validate_vcf(baseline_vcf, "Baseline")
    candidate_vcf = _validate_vcf(candidate_vcf, "Candidate")
    reference_sdf = os.path.abspath(os.fspath(reference_sdf))
    if not os.path.isdir(reference_sdf):
        raise FileNotFoundError(f"Reference SDF directory not found: {reference_sdf}")
    if regions_bed is not None:
        regions_bed = require_file(regions_bed, "Evaluation regions BED")

    output_dir = os.path.abspath(os.fspath(output_dir))
    baseline_output = os.path.join(output_dir, "baseline")
    candidate_output = os.path.join(output_dir, "candidate")
    report_file = os.path.join(output_dir, "benchmark.json")
    for path in (baseline_output, candidate_output, report_file):
        if os.path.exists(path):
            raise FileExistsError(f"Benchmark output already exists: {path}")

    executable = require_executable("rtg")
    os.makedirs(output_dir, exist_ok=True)
    baseline_metrics, baseline_summary = _run_vcfeval(
        executable,
        truth_vcf,
        baseline_vcf,
        reference_sdf,
        baseline_output,
        regions_bed,
    )
    candidate_metrics, candidate_summary = _run_vcfeval(
        executable,
        truth_vcf,
        candidate_vcf,
        reference_sdf,
        candidate_output,
        regions_bed,
    )
    delta = {
        metric: candidate_metrics[metric] - baseline_metrics[metric]
        for metric in (
            "true_positives_baseline",
            "true_positives_candidate",
            "false_positives",
            "false_negatives",
            "precision",
            "sensitivity",
            "f1_score",
        )
    }
    report = {
        "tool": "RTG vcfeval",
        "truth_vcf": truth_vcf,
        "reference_sdf": reference_sdf,
        "evaluation_regions_bed": regions_bed,
        "baseline": {
            "calls_vcf": baseline_vcf,
            "metrics": baseline_metrics,
            "summary_file": baseline_summary,
        },
        "candidate": {
            "calls_vcf": candidate_vcf,
            "metrics": candidate_metrics,
            "summary_file": candidate_summary,
        },
        "delta_candidate_minus_baseline": delta,
        "interpretation": (
            "Positive precision, sensitivity (recall), or F1 deltas indicate improvement "
            "relative to the supplied baseline call set."
        ),
        "report_file": report_file,
    }
    with open(report_file, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=2, sort_keys=True)
        file_handle.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Compare baseline and candidate VCFs against a truth VCF."
    )
    parser.add_argument("--truth", required=True, help="Trusted truth-set VCF/VCF.GZ")
    parser.add_argument("--baseline", required=True, help="Baseline caller VCF/VCF.GZ")
    parser.add_argument("--candidate", required=True, help="Candidate caller VCF/VCF.GZ")
    parser.add_argument(
        "--reference",
        required=True,
        help="RTG SDF reference matching the assembly used for both call sets",
    )
    parser.add_argument(
        "--regions",
        default=None,
        help="Optional BED of benchmarkable/confident regions",
    )
    parser.add_argument(
        "--output-dir",
        default="results/benchmark",
        help="Output directory (must not contain previous benchmark outputs)",
    )
    args = parser.parse_args()

    try:
        report = benchmark_variant_calls(
            args.truth,
            args.baseline,
            args.candidate,
            args.reference,
            args.output_dir,
            regions_bed=args.regions,
        )
    except (FileNotFoundError, FileExistsError, OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))

    for label in ("baseline", "candidate"):
        metrics = report[label]["metrics"]
        print(
            f"{label.capitalize()}: precision={metrics['precision']:.6f}, "
            f"sensitivity={metrics['sensitivity']:.6f}, "
            f"F1={metrics['f1_score']:.6f}"
        )
    print(f"Candidate minus baseline: {report['delta_candidate_minus_baseline']}")
    print(f"JSON report: {report['report_file']}")


if __name__ == "__main__":
    main()
