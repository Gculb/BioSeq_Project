import os
import re
import tempfile

from ._external import (
    default_output_path,
    require_executable,
    require_file,
    run_command,
    staged_output_path,
)


FLAGSTAT_LINE = re.compile(
    r"^\s*(\d+)\s+\+\s+(\d+)\s+(.+?)(?:\s+\(([^()]*)\))?\s*$"
)


def _resolve_samtools():
    return require_executable("samtools")


def quickcheck(alignment_file):
    """Check a SAM/BAM/CRAM file's header and EOF marker with samtools."""
    alignment_file = require_file(alignment_file, "Alignment")
    samtools = _resolve_samtools()
    run_command([samtools, "quickcheck", "-v", alignment_file])
    return True


def sort_bam(input_bam, output_bam=None, by_name=False):
    """Sort a BAM by coordinate or read name and return the output path."""
    input_bam = require_file(input_bam, "BAM")
    if output_bam is None:
        suffix = ".name_sorted.bam" if by_name else ".sorted.bam"
        output_bam = default_output_path(input_bam, suffix)
    output_bam = os.path.abspath(os.fspath(output_bam))
    if output_bam == input_bam:
        raise ValueError("output_bam must not overwrite the input BAM file.")

    samtools = _resolve_samtools()
    command = [samtools, "sort"]
    if by_name:
        command.append("-n")
    command.extend(["-O", "BAM", "-o"])
    with staged_output_path(output_bam) as temporary_path:
        command.extend([temporary_path, input_bam])
        run_command(command)
    return output_bam


def fixmate(input_bam, output_bam):
    """Add mate tags required by samtools markdup."""
    input_bam = require_file(input_bam, "Name-sorted BAM")
    output_bam = os.path.abspath(os.fspath(output_bam))
    if output_bam == input_bam:
        raise ValueError("output_bam must not overwrite the input BAM file.")
    samtools = _resolve_samtools()
    with staged_output_path(output_bam) as temporary_path:
        run_command(
            [samtools, "fixmate", "-m", "-O", "BAM", input_bam, temporary_path]
        )
    return output_bam


def mark_duplicates(input_bam, output_bam):
    """Mark duplicate reads in a coordinate-sorted BAM."""
    input_bam = require_file(input_bam, "Coordinate-sorted BAM")
    output_bam = os.path.abspath(os.fspath(output_bam))
    if output_bam == input_bam:
        raise ValueError("output_bam must not overwrite the input BAM file.")
    samtools = _resolve_samtools()
    with staged_output_path(output_bam) as temporary_path:
        run_command([samtools, "markdup", input_bam, temporary_path])
    return output_bam


def index_bam(bam_file):
    """Create a BAM index and return the absolute index path."""
    bam_file = require_file(bam_file, "BAM")
    index_file = f"{bam_file}.bai"
    samtools = _resolve_samtools()

    with staged_output_path(index_file) as temporary_path:
        run_command([samtools, "index", "-o", temporary_path, bam_file])
    return index_file


def _parse_percentage_pair(text):
    if not text or "%" not in text:
        return None
    percentages = []
    for item in text.split(":"):
        item = item.strip()
        match = re.fullmatch(r"(\d+(?:\.\d+)?)%", item)
        percentages.append(float(match.group(1)) if match else None)
    return percentages


def flagstat(bam_file):
    """Run samtools flagstat and return counts by category."""
    bam_file = require_file(bam_file, "BAM")
    samtools = _resolve_samtools()
    output = run_command([samtools, "flagstat", bam_file])
    metrics = {}

    for line in output.splitlines():
        match = FLAGSTAT_LINE.match(line)
        if match is None:
            continue
        passed, failed, category, annotation = match.groups()
        metrics[category] = {
            "passed": int(passed),
            "failed": int(failed),
        }
        percentages = _parse_percentage_pair(annotation)
        if percentages is not None:
            metrics[category]["percentages"] = percentages
        elif annotation:
            metrics[category]["description"] = annotation

    if not metrics:
        raise ValueError("samtools flagstat returned output that could not be parsed.")
    return metrics


def _parse_stat_value(value):
    value = value.strip()
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def stats(bam_file):
    """Run samtools stats and return parsed SN (summary) values."""
    bam_file = require_file(bam_file, "BAM")
    samtools = _resolve_samtools()
    output = run_command([samtools, "stats", bam_file])
    metrics = {}

    for line in output.splitlines():
        if not line.startswith("SN\t"):
            continue
        fields = line.split("\t", 2)
        if len(fields) < 3:
            continue
        name = fields[1].strip().rstrip(":")
        metrics[name] = _parse_stat_value(fields[2])

    if not metrics:
        raise ValueError("samtools stats returned no parseable summary (SN) records.")
    return metrics


def coverage(bam_file, region=None):
    """Return per-reference coverage metrics from ``samtools coverage``."""
    bam_file = require_file(bam_file, "BAM")
    samtools = _resolve_samtools()
    command = [samtools, "coverage"]
    if region is not None:
        if not isinstance(region, str) or not region.strip():
            raise ValueError("region must be a non-empty genomic interval.")
        command.extend(["-r", region])
    command.append(bam_file)
    output = run_command(command)

    header = None
    summaries = []
    metric_names = {
        "rname": "reference",
        "start": "start",
        "end": "end",
        "numreads": "read_count",
        "covbases": "covered_bases",
        "coverage": "coverage_percent",
        "meandepth": "mean_depth",
        "meanbaseq": "mean_base_quality",
        "meanmapq": "mean_mapping_quality",
    }
    integer_fields = {"start", "end", "numreads", "covbases"}
    float_fields = {"coverage", "meandepth", "meanbaseq", "meanmapq"}
    for line in output.splitlines():
        fields = line.split("\t")
        if fields[0].lstrip("#").casefold() == "rname":
            header = []
            for field in fields:
                field = field.lstrip("#").casefold()
                header.append({"startpos": "start", "endpos": "end"}.get(field, field))
            if header != list(metric_names):
                raise ValueError("samtools coverage returned an unexpected header.")
            continue
        if not fields or not fields[0].strip():
            continue
        if header is None or len(fields) != len(header):
            raise ValueError("samtools coverage returned an invalid summary row.")
        row = {}
        for key, value in zip(header, fields):
            if key == "rname":
                row[metric_names[key]] = value
            elif value == ".":
                row[metric_names[key]] = None
            elif key in integer_fields:
                row[metric_names[key]] = int(value)
            elif key in float_fields:
                row[metric_names[key]] = float(value)
        summaries.append(row)
    if not summaries:
        raise ValueError("samtools coverage returned no parseable summary rows.")
    return summaries


def depth(bam_file, output_file=None):
    """Run samtools depth, returning tab-separated text or writing it to a file.

    When ``output_file`` is supplied, output is written atomically and its path
    is returned. Otherwise the raw tab-separated depth table is returned.
    """
    bam_file = require_file(bam_file, "BAM")
    samtools = _resolve_samtools()
    command = [samtools, "depth", bam_file]
    if output_file is None:
        return run_command(command)

    output_file = os.path.abspath(os.fspath(output_file))
    if output_file == bam_file:
        raise ValueError("output_file must not overwrite the input BAM file.")
    return run_command(command, stdout_path=output_file)


def analyze_alignment(alignment_file, reference_file=None, region=None):
    """Validate and summarize a BAM or CRAM using samtools.

    CRAMs can use a local reference FASTA supplied with ``reference_file``.
    In that case the CRAM is temporarily decoded to BAM before its summaries
    are calculated, because samtools' summary commands do not consistently
    accept a per-command reference option.
    """
    alignment_file = require_file(alignment_file, "Alignment")
    extension = os.path.splitext(alignment_file)[1].lower()
    if extension not in {".bam", ".cram"}:
        raise ValueError("alignment_file must end in .bam or .cram.")

    quickcheck(alignment_file)
    analysis_file = alignment_file
    if extension == ".cram" and reference_file is not None:
        reference_file = require_file(reference_file, "Reference")
        samtools = _resolve_samtools()
        with tempfile.TemporaryDirectory(prefix="bioseq_cram_") as temp_dir:
            decoded_bam = os.path.join(temp_dir, "decoded.bam")
            with staged_output_path(decoded_bam) as temporary_path:
                run_command(
                    [
                        samtools,
                        "view",
                        "-T",
                        reference_file,
                        "-b",
                        "-o",
                        temporary_path,
                        alignment_file,
                    ]
                )
            metrics = {
                "flagstat": flagstat(decoded_bam),
                "stats": stats(decoded_bam),
                "coverage": coverage(decoded_bam, region=region),
            }
    else:
        metrics = {
            "flagstat": flagstat(analysis_file),
            "stats": stats(analysis_file),
            "coverage": coverage(analysis_file, region=region),
        }

    metrics["reference_file"] = (
        os.path.abspath(os.fspath(reference_file))
        if reference_file is not None
        else None
    )
    return metrics
