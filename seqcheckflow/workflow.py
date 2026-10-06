import hashlib
import json
import os
import re
import shutil
import subprocess
import threading
import time
from contextlib import contextmanager

import psutil

from ._external import require_executable, require_file, run_command
from .alignment import align_reads, convert_sam_to_bam, index_reference
from .qc import analyze_fastq_quality
from .samtools import (
    analyze_alignment,
    fixmate,
    index_bam,
    mark_duplicates,
    sort_bam,
)
from .validation import validate_sequence_file
from .variants import analyze_vcf


INTERVAL_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+:[0-9,]+-[0-9,]+$")


class _ProcessTreeSampler:
    def __init__(self, interval_seconds=0.05):
        self.root = psutil.Process()
        self.interval_seconds = interval_seconds
        self.stop_event = threading.Event()
        self.thread = None
        self.peak_rss_bytes = 0
        self.peak_process_count = 0
        self.cpu_seconds = 0.0
        self.start_cpu_seconds = 0.0
        self.process_cpu_by_id = {}

    def _sample(self):
        try:
            processes = [self.root, *self.root.children(recursive=True)]
        except psutil.Error:
            processes = [self.root]
        rss = 0
        cpu_seconds = 0.0
        process_count = 0
        for process in processes:
            try:
                rss += process.memory_info().rss
                cpu = process.cpu_times()
                process_id = (process.pid, process.create_time())
                process_cpu = cpu.user + cpu.system
                self.process_cpu_by_id[process_id] = max(
                    self.process_cpu_by_id.get(process_id, 0.0), process_cpu
                )
                process_count += 1
            except psutil.Error:
                continue
        cpu_seconds = sum(self.process_cpu_by_id.values())
        self.peak_rss_bytes = max(self.peak_rss_bytes, rss)
        self.peak_process_count = max(self.peak_process_count, process_count)
        self.cpu_seconds = cpu_seconds

    def _sample_until_stopped(self):
        while not self.stop_event.wait(self.interval_seconds):
            self._sample()


@contextmanager
def _monitor_process_tree():
    sampler = _ProcessTreeSampler()
    sampler._sample()
    sampler.start_cpu_seconds = sampler.cpu_seconds
    sampler.thread = threading.Thread(
        target=sampler._sample_until_stopped, daemon=True
    )
    sampler.thread.start()
    try:
        yield sampler
    finally:
        sampler.stop_event.set()
        sampler.thread.join()
        sampler._sample()
        sampler.process_tree_cpu_seconds = max(
            0.0, sampler.cpu_seconds - sampler.start_cpu_seconds
        )


def _profile_stage(stages, name, operation, output_paths=()):
    started = time.perf_counter()
    with _monitor_process_tree() as sampler:
        result = operation()
    stage = {
        "name": name,
        "wall_seconds": time.perf_counter() - started,
        "peak_process_tree_rss_bytes": sampler.peak_rss_bytes,
        "peak_process_count": sampler.peak_process_count,
        "process_tree_cpu_seconds": sampler.process_tree_cpu_seconds,
        "output_bytes": sum(
            os.path.getsize(path) for path in output_paths if os.path.isfile(path)
        ),
    }
    stages.append(stage)
    return result


def _sha256(filename):
    digest = hashlib.sha256()
    with open(filename, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tool_versions():
    commands = {
        "bwa": ["bwa"],
        "samtools": ["samtools", "--version"],
        "gatk": ["gatk", "--version"],
        "rtg": ["rtg", "version"],
    }
    versions = {}
    for name, command in commands.items():
        executable = require_executable(command[0])
        result = subprocess.run(
            [executable, *command[1:]],
            check=False,
            capture_output=True,
            text=True,
        )
        output = (result.stdout + "\n" + result.stderr).strip()
        if result.returncode not in {0, 1} or not output:
            raise RuntimeError(f"Could not determine {name} version: {output}")
        lines = [line.strip() for line in output.splitlines() if line.strip()]
        if name == "bwa":
            version_line = next(
                (line for line in lines if line.casefold().startswith("version:")),
                lines[0],
            )
            versions[name] = version_line
        else:
            versions[name] = lines[0]
    return versions


def _validate_interval(interval):
    if not isinstance(interval, str) or not INTERVAL_PATTERN.fullmatch(interval):
        raise ValueError("interval must have the form contig:start-end.")
    contig, coordinates = interval.split(":", 1)
    start, end = (int(value.replace(",", "")) for value in coordinates.split("-"))
    if start < 1 or end < start:
        raise ValueError("interval coordinates must be positive and ordered.")
    return contig, start, end


def _validate_fastqs(read1, read2):
    read1 = require_file(read1, "Read 1 FASTQ")
    read2 = require_file(read2, "Read 2 FASTQ")
    for filename in (read1, read2):
        validation = validate_sequence_file(filename, file_format="fastq")
        if not validation["valid"]:
            details = "\n".join(f"- {error}" for error in validation["errors"])
            raise ValueError(f"FASTQ validation failed for {filename}:\n{details}")
    return read1, read2


def _create_reference_copy(reference_fasta, work_dir):
    reference_copy = os.path.join(work_dir, "reference.fa")
    shutil.copy2(reference_fasta, reference_copy)
    return reference_copy


def _run_gatk_haplotype_caller(reference, alignment, interval, threads, output_vcf):
    gatk = require_executable("gatk")
    run_command(
        [
            gatk,
            "HaplotypeCaller",
            "-R",
            reference,
            "-I",
            alignment,
            "-O",
            output_vcf,
            "-L",
            interval,
            "--native-pair-hmm-threads",
            str(threads),
        ]
    )
    if not os.path.isfile(output_vcf):
        raise RuntimeError(f"GATK did not create its VCF output: {output_vcf}")
    return output_vcf


def run_germline_workflow(
    read1,
    read2,
    reference_fasta,
    known_sites_vcf,
    interval,
    output_dir,
    sample_name="HG002",
    threads=2,
    implementation="bioseq",
):
    """Run paired FASTQ alignment, duplicate marking, calling, and per-stage metrics.

    ``implementation="baseline"`` runs the documented BWA/samtools/GATK
    commands directly instead of using BioSeq's wrappers, for parity testing.
    """
    if implementation not in {"bioseq", "baseline"}:
        raise ValueError("implementation must be 'bioseq' or 'baseline'.")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer.")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample_name):
        raise ValueError("sample_name may contain only letters, numbers, _, ., and -.")
    _validate_interval(interval)

    read1, read2 = _validate_fastqs(read1, read2)
    reference_fasta = require_file(reference_fasta, "Reference FASTA")
    known_sites_vcf = require_file(known_sites_vcf, "Known-sites VCF")
    output_dir = os.path.abspath(os.fspath(output_dir))
    if os.path.exists(output_dir):
        raise FileExistsError(f"Workflow output directory already exists: {output_dir}")

    total_started = time.perf_counter()
    os.makedirs(output_dir)
    stages = []
    reference = _create_reference_copy(reference_fasta, output_dir)
    reference_index_files = [
        f"{reference}.{suffix}" for suffix in ("amb", "ann", "bwt", "pac", "sa")
    ]
    if implementation == "bioseq":
        _profile_stage(
            stages,
            "reference_index",
            lambda: index_reference(reference),
            reference_index_files,
        )
    else:
        bwa = require_executable("bwa")
        _profile_stage(
            stages,
            "reference_index",
            lambda: run_command([bwa, "index", reference]),
            reference_index_files,
        )

    reference_fai = f"{reference}.fai"
    reference_dict = os.path.splitext(reference)[0] + ".dict"
    if implementation == "bioseq":
        samtools = require_executable("samtools")
        gatk = require_executable("gatk")
        _profile_stage(
            stages,
            "reference_metadata",
            lambda: (
                run_command([samtools, "faidx", reference]),
                run_command(
                    [
                        gatk,
                        "CreateSequenceDictionary",
                        "-R",
                        reference,
                        "-O",
                        reference_dict,
                    ]
                ),
            ),
            [reference_fai, reference_dict],
        )
    else:
        samtools = require_executable("samtools")
        gatk = require_executable("gatk")

        def create_reference_metadata():
            run_command([samtools, "faidx", reference])
            run_command(
                [
                    gatk,
                    "CreateSequenceDictionary",
                    "-R",
                    reference,
                    "-O",
                    reference_dict,
                ]
            )

        _profile_stage(
            stages,
            "reference_metadata",
            create_reference_metadata,
            [reference_fai, reference_dict],
        )

    fastq_metrics = _profile_stage(
        stages,
        "fastq_qc",
        lambda: {
            "read1": analyze_fastq_quality(read1),
            "read2": analyze_fastq_quality(read2),
        },
    )
    sam_path = os.path.join(output_dir, "alignment.sam")
    bam_path = os.path.join(output_dir, "alignment.bam")
    name_sorted_bam = os.path.join(output_dir, "alignment.name_sorted.bam")
    fixmate_bam = os.path.join(output_dir, "alignment.fixmate.bam")
    coordinate_sorted_bam = os.path.join(output_dir, "alignment.sorted.bam")
    deduplicated_bam = os.path.join(output_dir, "alignment.deduplicated.bam")
    bam_index = f"{deduplicated_bam}.bai"
    read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

    if implementation == "bioseq":
        _profile_stage(
            stages,
            "bwa_mem_alignment",
            lambda: align_reads(
                [read1, read2],
                reference,
                sam_path,
                aligner="bwa",
                threads=threads,
                read_group=read_group,
            ),
            [sam_path],
        )
        _profile_stage(
            stages,
            "sam_to_bam",
            lambda: convert_sam_to_bam(sam_path, bam_path),
            [bam_path],
        )
        _profile_stage(
            stages,
            "name_sort",
            lambda: sort_bam(bam_path, name_sorted_bam, by_name=True),
            [name_sorted_bam],
        )
        _profile_stage(
            stages,
            "fixmate",
            lambda: fixmate(name_sorted_bam, fixmate_bam),
            [fixmate_bam],
        )
        _profile_stage(
            stages,
            "coordinate_sort",
            lambda: sort_bam(fixmate_bam, coordinate_sorted_bam),
            [coordinate_sorted_bam],
        )
        _profile_stage(
            stages,
            "duplicate_marking",
            lambda: mark_duplicates(coordinate_sorted_bam, deduplicated_bam),
            [deduplicated_bam],
        )
        _profile_stage(
            stages,
            "alignment_index",
            lambda: index_bam(deduplicated_bam),
            [bam_index],
        )
    else:
        bwa = require_executable("bwa")
        samtools = require_executable("samtools")

        def align():
            run_command(
                [
                    bwa,
                    "mem",
                    "-R",
                    read_group,
                    "-t",
                    str(threads),
                    reference,
                    read1,
                    read2,
                ],
                stdout_path=sam_path,
            )

        _profile_stage(stages, "bwa_mem_alignment", align, [sam_path])

        def convert():
            run_command([samtools, "view", "-b", "-o", bam_path, sam_path])

        _profile_stage(stages, "sam_to_bam", convert, [bam_path])

        def name_sort():
            run_command(
                [samtools, "sort", "-n", "-O", "BAM", "-o", name_sorted_bam, bam_path]
            )

        _profile_stage(stages, "name_sort", name_sort, [name_sorted_bam])

        def run_fixmate():
            run_command(
                [
                    samtools,
                    "fixmate",
                    "-m",
                    "-O",
                    "BAM",
                    name_sorted_bam,
                    fixmate_bam,
                ]
            )

        _profile_stage(stages, "fixmate", run_fixmate, [fixmate_bam])

        def coordinate_sort():
            run_command(
                [
                    samtools,
                    "sort",
                    "-O",
                    "BAM",
                    "-o",
                    coordinate_sorted_bam,
                    fixmate_bam,
                ]
            )

        _profile_stage(stages, "coordinate_sort", coordinate_sort, [coordinate_sorted_bam])

        def run_mark_duplicates():
            run_command(
                [samtools, "markdup", coordinate_sorted_bam, deduplicated_bam]
            )

        _profile_stage(
            stages, "duplicate_marking", run_mark_duplicates, [deduplicated_bam]
        )

        def create_alignment_index():
            run_command([samtools, "index", "-o", bam_index, deduplicated_bam])

        _profile_stage(stages, "alignment_index", create_alignment_index, [bam_index])

    alignment_metrics = _profile_stage(
        stages,
        "alignment_qc",
        lambda: analyze_alignment(deduplicated_bam, region=interval),
    )
    recalibration_table = os.path.join(output_dir, "bqsr.table")
    recalibrated_bam = os.path.join(output_dir, "alignment.recalibrated.bam")
    recalibrated_bam_index = f"{recalibrated_bam}.bai"
    gatk = require_executable("gatk")

    def run_base_recalibrator():
        run_command(
            [
                gatk,
                "BaseRecalibrator",
                "-R",
                reference,
                "-I",
                deduplicated_bam,
                "--known-sites",
                known_sites_vcf,
                "-L",
                interval,
                "-O",
                recalibration_table,
            ]
        )

    _profile_stage(
        stages,
        "base_recalibrator",
        run_base_recalibrator,
        [recalibration_table],
    )

    def run_apply_bqsr():
        run_command(
            [
                gatk,
                "ApplyBQSR",
                "-R",
                reference,
                "-I",
                deduplicated_bam,
                "--bqsr-recal-file",
                recalibration_table,
                "-L",
                interval,
                "-O",
                recalibrated_bam,
            ]
        )

    _profile_stage(stages, "apply_bqsr", run_apply_bqsr, [recalibrated_bam])
    if implementation == "bioseq":
        _profile_stage(
            stages,
            "recalibrated_alignment_index",
            lambda: index_bam(recalibrated_bam),
            [recalibrated_bam_index],
        )
    else:
        def index_recalibrated_alignment():
            run_command(
                [
                    require_executable("samtools"),
                    "index",
                    "-o",
                    recalibrated_bam_index,
                    recalibrated_bam,
                ]
            )

        _profile_stage(
            stages,
            "recalibrated_alignment_index",
            index_recalibrated_alignment,
            [recalibrated_bam_index],
        )
    output_vcf = os.path.join(output_dir, "calls.vcf.gz")
    _profile_stage(
        stages,
        "gatk_haplotype_caller",
        lambda: _run_gatk_haplotype_caller(
            reference, recalibrated_bam, interval, threads, output_vcf
        ),
        [output_vcf],
    )
    variant_result = _profile_stage(
        stages,
        "variant_qc",
        lambda: analyze_vcf(output_vcf),
    )
    if not variant_result["validation"]["valid"]:
        details = "\n".join(
            f"- {error}" for error in variant_result["validation"]["errors"]
        )
        raise ValueError(f"GATK produced an invalid VCF:\n{details}")
    tool_versions = _tool_versions()
    workflow_report = {
        "implementation": implementation,
        "sample_name": sample_name,
        "interval": interval,
        "threads": threads,
        "inputs": {
            "read1": {
                "path": read1,
                "size_bytes": os.path.getsize(read1),
                "sha256": _sha256(read1),
            },
            "read2": {
                "path": read2,
                "size_bytes": os.path.getsize(read2),
                "sha256": _sha256(read2),
            },
            "reference": {
                "path": reference_fasta,
                "size_bytes": os.path.getsize(reference_fasta),
                "sha256": _sha256(reference_fasta),
            },
            "known_sites_vcf": {
                "path": known_sites_vcf,
                "size_bytes": os.path.getsize(known_sites_vcf),
                "sha256": _sha256(known_sites_vcf),
            },
        },
        "tool_versions": tool_versions,
        "fastq_metrics": fastq_metrics,
        "alignment_metrics": alignment_metrics,
        "variant_validation": variant_result["validation"],
        "variant_metrics": variant_result["metrics"],
        "stages": stages,
        "total_wall_seconds": time.perf_counter() - total_started,
        "peak_stage_rss_bytes": max(
            (stage["peak_process_tree_rss_bytes"] for stage in stages), default=0
        ),
        "total_output_bytes": sum(stage["output_bytes"] for stage in stages),
        "reference_used": reference,
        "alignment_bam": deduplicated_bam,
        "recalibrated_alignment_bam": recalibrated_bam,
        "variant_vcf": output_vcf,
        "report_file": os.path.join(output_dir, "workflow.json"),
    }
    with open(workflow_report["report_file"], "w", encoding="utf-8") as file_handle:
        json.dump(workflow_report, file_handle, indent=2, sort_keys=True)
        file_handle.write("\n")
    return workflow_report
