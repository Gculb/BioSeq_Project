import argparse
import hashlib
import json
import os

from .benchmark import benchmark_variant_calls
from ._external import require_file
from .workflow import run_germline_workflow


PERFORMANCE_METRICS = (
    "wall_seconds",
    "peak_process_tree_rss_bytes",
    "process_tree_cpu_seconds",
    "output_bytes",
)


def _sha256_file(filename):
    digest = hashlib.sha256()
    with open(filename, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_directory(directory):
    digest = hashlib.sha256()
    found_file = False
    for root, directories, filenames in os.walk(directory):
        directories.sort()
        for filename in sorted(filenames):
            found_file = True
            path = os.path.join(root, filename)
            relative_path = os.path.relpath(path, directory).replace(os.sep, "/")
            digest.update(relative_path.encode("utf-8"))
            digest.update(_sha256_file(path).encode("ascii"))
    if not found_file:
        raise ValueError(f"Reference SDF directory contains no files: {directory}")
    return digest.hexdigest()


def compare_germline_workflows(
    read1,
    read2,
    reference_fasta,
    known_sites_vcf,
    reference_sdf,
    interval,
    truth_vcf,
    confident_regions_bed,
    output_dir,
    sample_name="HG002",
    threads=2,
):
    """Run direct-tool baseline and BioSeq wrapper workflows, then compare accuracy."""
    read1 = require_file(read1, "Read 1 FASTQ")
    read2 = require_file(read2, "Read 2 FASTQ")
    reference_fasta = require_file(reference_fasta, "Reference FASTA")
    known_sites_vcf = require_file(known_sites_vcf, "Known-sites VCF")
    truth_vcf = require_file(truth_vcf, "Truth VCF")
    confident_regions_bed = require_file(
        confident_regions_bed, "Confident-regions BED"
    )
    reference_sdf = os.path.abspath(os.fspath(reference_sdf))
    if not os.path.isdir(reference_sdf):
        raise FileNotFoundError(f"Reference SDF directory not found: {reference_sdf}")
    output_dir = os.path.abspath(os.fspath(output_dir))
    if os.path.exists(output_dir):
        raise FileExistsError(f"Benchmark output directory already exists: {output_dir}")
    os.makedirs(output_dir)

    baseline = run_germline_workflow(
        read1,
        read2,
        reference_fasta,
        known_sites_vcf,
        interval,
        os.path.join(output_dir, "baseline"),
        sample_name=sample_name,
        threads=threads,
        implementation="baseline",
    )
    bioseq = run_germline_workflow(
        read1,
        read2,
        reference_fasta,
        known_sites_vcf,
        interval,
        os.path.join(output_dir, "bioseq"),
        sample_name=sample_name,
        threads=threads,
        implementation="bioseq",
    )
    if baseline["tool_versions"] != bioseq["tool_versions"]:
        raise RuntimeError("Baseline and BioSeq runs used different tool versions.")

    accuracy = benchmark_variant_calls(
        truth_vcf,
        baseline["variant_vcf"],
        bioseq["variant_vcf"],
        reference_sdf,
        os.path.join(output_dir, "truth_evaluation"),
        regions_bed=confident_regions_bed,
    )
    baseline_stages = {stage["name"]: stage for stage in baseline["stages"]}
    bioseq_stages = {stage["name"]: stage for stage in bioseq["stages"]}
    if baseline_stages.keys() != bioseq_stages.keys():
        raise RuntimeError("Baseline and BioSeq runs executed different stages.")

    stage_deltas = {}
    for name in baseline_stages:
        stage_deltas[name] = {
            metric: bioseq_stages[name][metric] - baseline_stages[name][metric]
            for metric in PERFORMANCE_METRICS
        }
    report_path = os.path.join(output_dir, "comparison.json")
    report = {
        "sample": sample_name,
        "region": interval,
        "threads": threads,
        "tool_versions": bioseq["tool_versions"],
        "benchmark_inputs": {
            "truth_vcf": {
                "path": truth_vcf,
                "sha256": _sha256_file(truth_vcf),
            },
            "confident_regions_bed": {
                "path": confident_regions_bed,
                "sha256": _sha256_file(confident_regions_bed),
            },
            "reference_sdf": {
                "path": reference_sdf,
                "sha256": _sha256_directory(reference_sdf),
            },
        },
        "baseline": {
            "definition": (
                "Direct BWA, samtools, and GATK CLI calls using the same fixed "
                "arguments as the BioSeq wrapper workflow."
            ),
            "workflow_report": baseline["report_file"],
            "total_wall_seconds": baseline["total_wall_seconds"],
            "peak_stage_rss_bytes": baseline["peak_stage_rss_bytes"],
            "total_output_bytes": baseline["total_output_bytes"],
            "stages": baseline_stages,
        },
        "bioseq": {
            "definition": "Same tool versions and arguments invoked through BioSeq wrappers.",
            "workflow_report": bioseq["report_file"],
            "total_wall_seconds": bioseq["total_wall_seconds"],
            "peak_stage_rss_bytes": bioseq["peak_stage_rss_bytes"],
            "total_output_bytes": bioseq["total_output_bytes"],
            "stages": bioseq_stages,
        },
        "candidate_minus_baseline": {
            "total_wall_seconds": (
                bioseq["total_wall_seconds"] - baseline["total_wall_seconds"]
            ),
            "peak_stage_rss_bytes": (
                bioseq["peak_stage_rss_bytes"] - baseline["peak_stage_rss_bytes"]
            ),
            "total_output_bytes": (
                bioseq["total_output_bytes"] - baseline["total_output_bytes"]
            ),
            "stages": stage_deltas,
        },
        "truth_evaluation": accuracy,
        "measurement_notes": [
            "Each workflow is run once; rerun in a fresh output directory to estimate timing variability.",
            "Peak process-tree RSS is sampled every 50 ms and may miss short-lived peaks.",
            "Process-tree RSS is summed across processes and can double-count shared memory.",
            "The timer includes Python and subprocess execution but excludes data preparation and truth-set downloads.",
            "This regional pilot assesses small-variant calls only and is not clinical validation.",
        ],
        "report_file": report_path,
    }
    with open(report_path, "w", encoding="utf-8") as file_handle:
        json.dump(report, file_handle, indent=2, sort_keys=True)
        file_handle.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Run BioSeq and direct-tool germline workflows and compare to GIAB truth."
    )
    parser.add_argument("--read1", required=True)
    parser.add_argument("--read2", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--known-sites-vcf", required=True)
    parser.add_argument("--reference-sdf", required=True)
    parser.add_argument("--truth-vcf", required=True)
    parser.add_argument("--confident-regions", required=True)
    parser.add_argument("--interval", default="chr20:10000000-11000000")
    parser.add_argument("--sample-name", default="HG002")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output-dir", default="results/HG002_chr20_comparison")
    args = parser.parse_args()

    try:
        report = compare_germline_workflows(
            args.read1,
            args.read2,
            args.reference_fasta,
            args.known_sites_vcf,
            args.reference_sdf,
            args.interval,
            args.truth_vcf,
            args.confident_regions,
            args.output_dir,
            sample_name=args.sample_name,
            threads=args.threads,
        )
    except (FileNotFoundError, FileExistsError, OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))

    comparison = report["candidate_minus_baseline"]
    print(
        "BioSeq minus direct-tool baseline: "
        f"{comparison['total_wall_seconds']:.3f} wall seconds, "
        f"{comparison['peak_stage_rss_bytes']} peak-RSS bytes"
    )
    print(f"Comparison report: {report['report_file']}")


if __name__ == "__main__":
    main()
