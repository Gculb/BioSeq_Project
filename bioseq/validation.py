import os

from .fastq import read_fastq
from .variants import analyze_vcf


NUCLEOTIDE_ALPHABET = set("ACGTUNRYKMSWBDHV")
FASTQ_EXTENSIONS = {".fastq", ".fq"}
FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}
ALIGNMENT_EXTENSIONS = {".bam", ".cram"}


def detect_sequence_format(filename):
    """Return the format indicated by a supported sequence-file extension."""
    filename = os.fspath(filename).lower()
    basename = os.path.basename(filename)
    if basename.endswith(".vcf.gz"):
        return "vcf"
    extension = os.path.splitext(filename)[1]
    if extension in FASTA_EXTENSIONS:
        return "fasta"
    if extension in FASTQ_EXTENSIONS:
        return "fastq"
    if extension == ".vcf":
        return "vcf"
    if extension in ALIGNMENT_EXTENSIONS:
        return extension[1:]
    raise ValueError(f"Unsupported sequence file extension: {extension or '(none)'}")


def validate_sequence_file(filename, file_format=None):
    """Validate sequence records and return a structured validation summary."""
    filename = os.fspath(filename)
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Sequence file not found: {filename}")

    if file_format is None:
        file_format = detect_sequence_format(filename)
    elif file_format not in {"fasta", "fastq", "vcf", "bam", "cram"}:
        raise ValueError("file_format must be fasta, fastq, vcf, bam, or cram.")

    errors = []
    record_count = 0
    try:
        if file_format == "fastq":
            records = read_fastq(filename)
            record_count = len(records)
            if not records:
                errors.append("File contains no FASTQ records.")
            for index, record in enumerate(records, start=1):
                invalid_bases = sorted(set(record["sequence"].upper()) - NUCLEOTIDE_ALPHABET)
                if not record["sequence"]:
                    errors.append(f"FASTQ record {index} has an empty sequence.")
                if invalid_bases:
                    errors.append(
                        f"FASTQ record {index} contains invalid nucleotide characters: "
                        f"{''.join(invalid_bases)}."
                    )
                invalid_quality = [
                    character
                    for character in record["quality"]
                    if not 33 <= ord(character) <= 126
                ]
                if invalid_quality:
                    errors.append(
                        f"FASTQ record {index} contains characters outside the "
                        "Phred+33 quality range."
                    )
        elif file_format == "fasta":
            record_count, fasta_errors = _validate_fasta(filename)
            errors.extend(fasta_errors)
        elif file_format == "vcf":
            result = analyze_vcf(filename)
            return result["validation"]
        else:
            from .samtools import quickcheck

            quickcheck(filename)
            record_count = None
    except (UnicodeError, ValueError) as exc:
        errors.append(str(exc))
    except RuntimeError as exc:
        errors.append(str(exc))

    if record_count == 0 and not errors and file_format in {"fasta", "fastq"}:
        errors.append(f"File contains no {file_format.upper()} records.")

    return {
        "valid": not errors,
        "format": file_format,
        "record_count": record_count,
        "errors": errors,
    }


def _validate_fasta(filename):
    errors = []
    record_count = 0
    current_id = None
    current_length = 0
    identifiers = set()

    def finish_record(line_number):
        if current_id is not None and current_length == 0:
            errors.append(f"FASTA record '{current_id}' has no sequence before line {line_number}.")

    with open(filename, "r", encoding="utf-8") as file_handle:
        for line_number, raw_line in enumerate(file_handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                finish_record(line_number)
                header = line[1:].strip()
                if not header:
                    errors.append(f"FASTA header on line {line_number} is empty.")
                    current_id = None
                    current_length = 0
                    continue

                current_id = header.split()[0]
                current_length = 0
                record_count += 1
                if current_id in identifiers:
                    errors.append(f"Duplicate FASTA identifier '{current_id}'.")
                identifiers.add(current_id)
                continue

            if current_id is None:
                errors.append(f"Sequence data appears before a FASTA header on line {line_number}.")
                continue

            sequence = line.upper()
            invalid_bases = sorted(set(sequence) - NUCLEOTIDE_ALPHABET)
            if invalid_bases:
                errors.append(
                    f"Invalid nucleotide characters on FASTA line {line_number}: "
                    f"{''.join(invalid_bases)}."
                )
            current_length += len(sequence)

    finish_record("end of file")
    if record_count == 0:
        errors.append("File contains no FASTA records.")

    return record_count, errors
