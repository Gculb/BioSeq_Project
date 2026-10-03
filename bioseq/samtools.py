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


def sort_bam(input_bam, output_bam=None):
    """Sort a BAM with samtools and return the absolute sorted BAM path."""
    input_bam = require_file(input_bam, "BAM")
    if output_bam is None:
        output_bam = default_output_path(input_bam, ".sorted.bam")
    output_bam = os.path.abspath(os.fspath(output_bam))
    if output_bam == input_bam:
        raise ValueError("output_bam must not overwrite the input BAM file.")

    samtools = _resolve_samtools()
    with staged_output_path(output_bam) as temporary_path:
        run_command(
            [samtools, "sort", "-O", "BAM", "-o", temporary_path, input_bam]
        )
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


def analyze_alignment(alignment_file, reference_file=None):
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
            }
    else:
        metrics = {
            "flagstat": flagstat(analysis_file),
            "stats": stats(analysis_file),
        }

    metrics["reference_file"] = (
        os.path.abspath(os.fspath(reference_file))
        if reference_file is not None
        else None
    )
    return metrics
