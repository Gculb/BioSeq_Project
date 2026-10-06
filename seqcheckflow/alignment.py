import os

from ._external import (
    default_output_path,
    require_executable,
    require_file,
    run_command,
    staged_output_path,
)


def align_reads(
    reads_file,
    reference_file,
    output_sam,
    aligner="bwa",
    threads=1,
    read_group=None,
):
    """Align one or two FASTQ files to a reference with BWA-MEM or minimap2.

    Returns the absolute path to the generated SAM file. BWA references must
    first be indexed with ``bwa index reference.fasta``; minimap2 builds or
    consumes its index automatically. Pass two files for paired-end reads.
    """
    if isinstance(reads_file, (str, os.PathLike)):
        reads_files = [reads_file]
    elif reads_file is None:
        raise ValueError("reads_file must be one FASTQ path or a pair of FASTQ paths.")
    else:
        reads_files = list(reads_file)
    if len(reads_files) not in {1, 2}:
        raise ValueError("reads_file must be one FASTQ path or a pair of FASTQ paths.")

    reads_files = [require_file(path, "Reads") for path in reads_files]
    reference_file = require_file(reference_file, "Reference")
    output_sam = os.path.abspath(os.fspath(output_sam))
    if output_sam in {*reads_files, reference_file}:
        raise ValueError("output_sam must not overwrite the reads or reference file.")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer.")

    if aligner == "bwa":
        executable = require_executable("bwa")
        if not os.path.isfile(f"{reference_file}.bwt"):
            raise FileNotFoundError(
                f"BWA index not found for {reference_file}; create it with "
                f"`bwa index {reference_file}`."
            )
        command = [
            executable,
            "mem",
            "-t",
            str(threads),
            reference_file,
            *reads_files,
        ]
        if read_group is not None:
            if not isinstance(read_group, str) or not read_group.startswith("@RG\\t"):
                raise ValueError("read_group must be an @RG header line.")
            command[2:2] = ["-R", read_group]
    elif aligner == "minimap2":
        executable = require_executable("minimap2")
        command = [
            executable,
            "-a",
            "-t",
            str(threads),
            reference_file,
            *reads_files,
        ]
    else:
        raise ValueError("aligner must be 'bwa' or 'minimap2'.")

    return run_command(command, stdout_path=output_sam)


def index_reference(reference_file):
    """Create a BWA index for a reference FASTA."""
    reference_file = require_file(reference_file, "Reference")
    executable = require_executable("bwa")
    run_command([executable, "index", reference_file])
    return reference_file


def convert_sam_to_bam(sam_file, output_bam=None):
    """Convert a SAM file to BAM with samtools and return the BAM path."""
    sam_file = require_file(sam_file, "SAM")
    if output_bam is None:
        output_bam = default_output_path(sam_file, ".bam")
    output_bam = os.path.abspath(os.fspath(output_bam))
    if output_bam == sam_file:
        raise ValueError("output_bam must not overwrite the input SAM file.")

    samtools = require_executable("samtools")
    command = [samtools, "view", "-b", "-o"]
    with staged_output_path(output_bam) as temporary_path:
        command.extend([temporary_path, sam_file])
        run_command(command)
    return output_bam
