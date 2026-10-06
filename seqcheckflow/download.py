import argparse
import glob
import os
import re
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.parse
import urllib.request

from Bio import Entrez


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DATA_DIR = os.path.join(PROJECT_ROOT, "data")
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
DOWNLOAD_TIMEOUT_SECONDS = 60


def _fastq_output_paths(output_dir, run_accession):
    pattern = re.compile(rf"{re.escape(run_accession)}(?:_\d+)?\.fastq$")
    candidates = glob.glob(os.path.join(output_dir, f"{run_accession}*.fastq"))
    return sorted(
        path for path in candidates if pattern.fullmatch(os.path.basename(path))
    )


def _resolve_output_dir(output_dir=None):
    if output_dir is None:
        output_dir = DEFAULT_DATA_DIR
    elif not os.path.isabs(output_dir):
        output_dir = os.path.join(PROJECT_ROOT, output_dir)

    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def _download_url_file(url, output_dir, filename, allowed_extensions, file_type):
    if not isinstance(url, str):
        raise ValueError(f"{file_type} URL must be a string.")
    url = url.strip()
    parsed_url = urllib.parse.urlsplit(url)
    if parsed_url.scheme not in {"http", "https"} or not parsed_url.hostname:
        raise ValueError(f"{file_type} URL must be an absolute HTTP or HTTPS URL.")
    if parsed_url.username or parsed_url.password:
        raise ValueError(f"{file_type} URL must not include embedded credentials.")

    if filename is None:
        filename = os.path.basename(urllib.parse.unquote(parsed_url.path))
    if not filename or filename in {".", ".."} or os.path.basename(filename) != filename:
        raise ValueError("filename must be a plain file name without directory paths.")

    normalized_filename = filename.lower()
    if not any(normalized_filename.endswith(ext) for ext in allowed_extensions):
        allowed = ", ".join(allowed_extensions)
        raise ValueError(f"{file_type} filename must end with one of: {allowed}.")

    output_dir = _resolve_output_dir(output_dir)
    output_path = os.path.abspath(os.path.join(output_dir, filename))
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file already exists: {output_path}")

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "SeqCheckFlow/1.0 (sequence data downloader)"},
    )
    descriptor, temporary_path = tempfile.mkstemp(
        prefix=f".{filename}.", suffix=".part", dir=output_dir
    )
    os.close(descriptor)
    downloaded_bytes = 0
    try:
        with urllib.request.urlopen(
            request, timeout=DOWNLOAD_TIMEOUT_SECONDS
        ) as response, open(temporary_path, "wb") as output:
            while chunk := response.read(DOWNLOAD_CHUNK_SIZE):
                output.write(chunk)
                downloaded_bytes += len(chunk)

        if not downloaded_bytes:
            raise ValueError(f"Downloaded {file_type} file is empty: {url}")

        _validate_downloaded_file(temporary_path, filename, file_type)
        if os.path.exists(output_path):
            raise FileExistsError(f"Output file already exists: {output_path}")
        os.replace(temporary_path, output_path)
    except urllib.error.URLError as exc:
        reason = getattr(exc, "reason", exc)
        raise RuntimeError(f"Failed to download {file_type} file from {url}: {reason}") from exc
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)

    return output_path


def _validate_downloaded_file(file_path, filename, file_type):
    lowered_name = filename.lower()
    if file_type == "VCF":
        if lowered_name.endswith(".vcf.gz"):
            import gzip

            opener = gzip.open
        else:
            opener = open

        try:
            with opener(file_path, "rb") as file_handle:
                header = file_handle.readline(1024)
        except (OSError, EOFError) as exc:
            raise ValueError(f"Downloaded file is not a readable VCF: {filename}") from exc

        if not header.startswith(b"##fileformat=VCF"):
            raise ValueError(f"Downloaded file does not have a VCF header: {filename}")
        return

    if lowered_name.endswith(".cram"):
        with open(file_path, "rb") as file_handle:
            signature = file_handle.read(4)
        if signature != b"CRAM":
            raise ValueError(f"Downloaded file does not have a CRAM signature: {filename}")
        return

    import gzip

    try:
        with gzip.open(file_path, "rb") as file_handle:
            signature = file_handle.read(4)
    except (OSError, EOFError) as exc:
        raise ValueError(f"Downloaded file is not a valid compressed BAM: {filename}") from exc
    if signature != b"BAM\x01":
        raise ValueError(f"Downloaded file does not have a BAM signature: {filename}")


def download_variant_vcf(url, output_dir=None, filename=None):
    """Download a VCF or compressed VCF from a direct HTTP(S) URL.

    By default, variant files are saved in ``data/variants/``.
    """
    return _download_url_file(
        url,
        output_dir or os.path.join(DEFAULT_DATA_DIR, "variants"),
        filename,
        (".vcf", ".vcf.gz"),
        "VCF",
    )


def download_variants(urls, output_dir=None):
    """Download multiple VCF files and return their paths."""
    if isinstance(urls, str):
        urls = [urls]
    if urls is None:
        return []

    return [
        download_variant_vcf(url, output_dir=output_dir)
        for url in urls
        if url is not None and str(url).strip()
    ]


def download_alignment_file(url, output_dir=None, filename=None):
    """Download an aligned BAM or CRAM file from a direct HTTP(S) URL.

    By default, alignments are saved in ``data/alignments/``. BAM files may
    be coordinate sorted and indexed afterward using ``seqcheckflow.samtools``.
    """
    return _download_url_file(
        url,
        output_dir or os.path.join(DEFAULT_DATA_DIR, "alignments"),
        filename,
        (".bam", ".cram"),
        "alignment",
    )


def download_bam(url, output_dir=None, filename=None):
    """Download a BAM file from a direct HTTP(S) URL."""
    return _download_url_file(
        url,
        output_dir or os.path.join(DEFAULT_DATA_DIR, "alignments"),
        filename,
        (".bam",),
        "BAM",
    )


def search_sequences(term, email, database="nucleotide", retmax=10):
    """Search NCBI for matching sequence records and return matching ID strings."""
    term = str(term).strip()
    if not term:
        raise ValueError("A search term is required.")

    if not email or not str(email).strip():
        raise ValueError("An email address is required for NCBI Entrez.")

    Entrez.email = email
    handle = Entrez.esearch(db=database, term=term, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()

    return result.get("IdList", [])


def download_sequence(accession, email, output_dir=None, database="nucleotide"):
    """Download one sequence from NCBI and save it to the data directory."""
    accession = str(accession).strip()
    if not accession:
        raise ValueError("Accession is required.")

    if not email or not str(email).strip():
        raise ValueError("An email address is required for NCBI Entrez.")

    output_dir = _resolve_output_dir(output_dir)
    Entrez.email = email

    handle = Entrez.efetch(db=database, id=accession, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()

    if not fasta_data.strip():
        raise ValueError(f"No sequence returned for accession {accession}.")

    file_path = os.path.join(output_dir, f"{accession}.fasta")
    with open(file_path, "w", encoding="utf-8") as file_handle:
        file_handle.write(fasta_data)

    return file_path


def download_sequences(accessions, email, output_dir=None, database="nucleotide"):
    """Download multiple sequences and return a list of saved file paths."""
    if accessions is None:
        return []

    if isinstance(accessions, str):
        accessions = [accessions]

    downloaded = []
    for accession in accessions:
        if accession is None or not str(accession).strip():
            continue
        downloaded.append(download_sequence(accession, email, output_dir=output_dir, database=database))

    return downloaded


def download_fastq(run_accession, output_dir=None, threads=1):
    """Download an SRA run as one or more FASTQ files using fasterq-dump.

    SRA Toolkit must be installed and ``fasterq-dump`` must be on PATH.
    Paired-end runs are written to separate ``_1`` and ``_2`` FASTQ files.
    """
    run_accession = str(run_accession).strip()
    if not re.fullmatch(r"(?:SRR|ERR|DRR)\d+", run_accession, re.IGNORECASE):
        raise ValueError("Provide a valid SRA run accession, such as SRR123456.")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer.")

    fasterq_dump = shutil.which("fasterq-dump")
    if fasterq_dump is None:
        raise FileNotFoundError(
            "fasterq-dump was not found. Install NCBI SRA Toolkit and add it to PATH."
        )

    output_dir = _resolve_output_dir(output_dir)
    if _fastq_output_paths(output_dir, run_accession):
        raise FileExistsError(
            f"FASTQ output for {run_accession} already exists in {output_dir}."
        )

    command = [
        fasterq_dump,
        "--split-files",
        "--threads",
        str(threads),
        "--outdir",
        output_dir,
        run_accession,
    ]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or "").strip()
        message = f"fasterq-dump failed for {run_accession}"
        if details:
            message += f": {details}"
        raise RuntimeError(message) from exc

    output_files = _fastq_output_paths(output_dir, run_accession)
    if not output_files:
        raise RuntimeError(f"fasterq-dump produced no FASTQ files for {run_accession}.")

    return output_files


def download_fastqs(run_accessions, output_dir=None, threads=1):
    """Download multiple SRA runs and return their FASTQ file paths."""
    if isinstance(run_accessions, str):
        run_accessions = [run_accessions]

    if run_accessions is None:
        return []

    downloaded = []
    for run_accession in run_accessions:
        if run_accession is None or not str(run_accession).strip():
            continue
        downloaded.extend(
            download_fastq(run_accession, output_dir=output_dir, threads=threads)
        )
    return downloaded


def main():
    parser = argparse.ArgumentParser(
        description="Search for and download sequence data from NCBI into the local data directory."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    search_parser = subparsers.add_parser("search", help="Find sequence IDs matching a search term")
    search_parser.add_argument("term", help="NCBI search term, for example 'Escherichia coli genome'")
    search_parser.add_argument("--email", required=True, help="Your email for NCBI Entrez requests")
    search_parser.add_argument("--db", default="nucleotide", help="NCBI database to search (default: nucleotide)")
    search_parser.add_argument("--retmax", type=int, default=10, help="Maximum number of results to return")

    download_parser = subparsers.add_parser("download", help="Download one or more accessions to the data directory")
    download_parser.add_argument("accessions", nargs="+", help="One or more NCBI accession numbers")
    download_parser.add_argument("--email", required=True, help="Your email for NCBI Entrez requests")
    download_parser.add_argument("--db", default="nucleotide", help="NCBI database to fetch (default: nucleotide)")
    download_parser.add_argument(
        "--dir",
        default="data",
        help="Destination directory relative to the project root (default: data)",
    )

    fastq_parser = subparsers.add_parser(
        "fastq", help="Download one or more SRA runs as FASTQ files"
    )
    fastq_parser.add_argument(
        "run_accessions", nargs="+", help="SRA run accession(s), such as SRR123456"
    )
    fastq_parser.add_argument(
        "--dir",
        default="data",
        help="Destination directory relative to the project root (default: data)",
    )
    fastq_parser.add_argument(
        "--threads", type=int, default=1, help="Number of fasterq-dump threads"
    )

    variants_parser = subparsers.add_parser(
        "variants", help="Download one or more VCF files from direct URLs"
    )
    variants_parser.add_argument(
        "urls", nargs="+", help="HTTP(S) URLs ending in .vcf or .vcf.gz"
    )
    variants_parser.add_argument(
        "--dir",
        default="data/variants",
        help="Destination directory relative to the project root",
    )

    bam_parser = subparsers.add_parser(
        "bam", help="Download one or more BAM files from direct URLs"
    )
    bam_parser.add_argument(
        "urls", nargs="+", help="HTTP(S) URLs ending in .bam"
    )
    bam_parser.add_argument(
        "--dir",
        default="data/alignments",
        help="Destination directory relative to the project root",
    )

    alignment_parser = subparsers.add_parser(
        "alignment", help="Download one or more BAM or CRAM files from direct URLs"
    )
    alignment_parser.add_argument(
        "urls", nargs="+", help="HTTP(S) URLs ending in .bam or .cram"
    )
    alignment_parser.add_argument(
        "--dir",
        default="data/alignments",
        help="Destination directory relative to the project root",
    )

    args = parser.parse_args()

    if args.command == "search":
        ids = search_sequences(args.term, args.email, database=args.db, retmax=args.retmax)
        if not ids:
            print("No matching records found.")
            return

        print("\n".join(ids))
        return

    if args.command == "download":
        paths = download_sequences(args.accessions, args.email, output_dir=args.dir, database=args.db)
        for path in paths:
            print(path)
        return

    if args.command == "fastq":
        paths = download_fastqs(
            args.run_accessions, output_dir=args.dir, threads=args.threads
        )
        for path in paths:
            print(path)
        return

    if args.command == "variants":
        paths = download_variants(args.urls, output_dir=args.dir)
        for path in paths:
            print(path)
        return

    if args.command == "bam":
        for url in args.urls:
            print(download_bam(url, output_dir=args.dir))
        return

    if args.command == "alignment":
        for url in args.urls:
            print(download_alignment_file(url, output_dir=args.dir))
        return

    parser.error("Unknown command")


if __name__ == "__main__":
    main()
