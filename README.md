# BioSeq Starter

A small Python bioinformatics pipeline for reasoning about sequence data.

## Current version

The project includes FASTA/FASTQ validation and analysis, FASTQ quality control,
NCBI FASTA and SRA FASTQ download helpers, and JSON/text pipeline reports.
The local `bioseq.pipeline` processes **one input file at a time**. It accepts
FASTA, FASTQ, VCF/VCF.GZ, BAM, and CRAM files and writes format-specific
quality-control/statistics reports. It does not automatically download data,
run BLAST, align reads, or call variants.

| File type | Can download? | Accepted by `bioseq.pipeline`? | Available processing |
| --- | --- | --- | --- |
| FASTA | Yes | Yes | Validate, sequence counts/lengths, GC content, JSON/text report |
| FASTQ | Yes | Yes | Validate, read-length and Phred quality metrics, GC/ambiguous-base metrics, JSON/text report |
| BAM | Yes, from a direct URL | No | Separate samtools wrappers can sort, index, and summarize BAMs |
| CRAM | Yes, from a direct URL | No | No CRAM-specific pipeline workflow; a reference may be needed by tools |
| VCF / VCF.GZ | Yes, from a direct URL | Yes | Variant-record, allele, filter, sample-genotype, and QUAL summaries; no annotation or interpretation |
| BAM | Yes, from a direct URL | Yes | `samtools quickcheck`, `flagstat`, and `stats` summaries |
| CRAM | Yes, from a direct URL | Yes | Same samtools summaries; may need a reference FASTA |

The `bioseq.pipeline` analysis command does not download files automatically.
Download and analysis remain separate steps.

Optional external-tool methods are available in `bioseq/blast.py`,
`bioseq/alignment.py`, and `bioseq/samtools.py`. `BLAST_search` submits a remote
query to NCBI; it requires an internet connection and is subject to NCBI usage
limits. Read alignment requires BWA (`bwa`) or minimap2 (`minimap2`), while
SAM/BAM conversion and BAM operations require samtools (`samtools`). These
executables must be installed separately and available on `PATH`.

The alignment wrapper is a separate step from `bioseq.pipeline`. For BWA,
index a reference once before aligning:

```bash
bwa index data/examples/example.fasta
```

The Python wrappers return result paths for file-producing methods and parsed
metrics for `flagstat`/`stats`. Paired-end reads can be passed to `align_reads`
as a two-item list or tuple.

Example Python usage:

```python
from bioseq.alignment import align_reads, convert_sam_to_bam
from bioseq.samtools import flagstat, index_bam, sort_bam

sam = align_reads(
    ["data/raw/sample_1.fastq", "data/raw/sample_2.fastq"],
    "data/references/reference.fasta",
    "data/processed/sample.sam",
    aligner="bwa",
    threads=4,
)
bam = convert_sam_to_bam(sam)
sorted_bam = sort_bam(bam)
index_bam(sorted_bam)
mapping_summary = flagstat(sorted_bam)
```

For a remote BLAST query:

```python
from bioseq.blast import BLAST_search

hits = BLAST_search(
    "ATGCGTACGT",
    program="blastn",
    database="nt",
    email="you@example.com",
)
```

## Run

From this directory:

```bash
python -m bioseq.cli data/examples/example.fasta
```

or:

```bash
python -m bioseq.cli data/examples/example.fastq
```

Run the local validation, analysis, and reporting workflow:

```bash
python -m bioseq.pipeline data/examples/example.fastq
```

Analyze other supported file types the same way:

```bash
python -m bioseq.pipeline data/variants/sample.vcf.gz
python -m bioseq.pipeline data/alignments/sample.bam
python -m bioseq.pipeline data/alignments/sample.cram --reference data/references/reference.fasta
```

Analyzing BAM/CRAM requires samtools installed on `PATH`. A supplied reference
allows CRAM to be decoded through a temporary BAM while summary metrics are
calculated. Without one, samtools must be able to resolve the reference from
the CRAM metadata or its configured reference cache.

The repository includes small synthetic FASTA and FASTQ examples under
`data/examples/`; input files in `data/raw/` remain ignored and local.

Reports are written to `results/` by default. Override the destination or the
mean-read Phred threshold used for FASTQ low-quality percentages with:

```bash
python -m bioseq.pipeline data/examples/example.fastq \
  --results-dir results \
  --minimum-mean-quality 20
```

The equivalent convenience script is `python scripts/run_pipeline.py <input>`.

## Download sequence data

Search NCBI and download nucleotide records as FASTA:

```bash
python -m bioseq.download search "Escherichia coli" --email you@example.com
python -m bioseq.download download NC_000913 --email you@example.com
```

Download SRA runs as FASTQ (requires the NCBI SRA Toolkit, including
`fasterq-dump`, installed and available on `PATH`):

```bash
python -m bioseq.download fastq SRR123456 --threads 4
```

Downloads are saved in `data/` by default. Paired-end runs produce separate
`_1.fastq` and `_2.fastq` files.

Download variant calls as VCF/VCF.GZ and alignment files as BAM or CRAM from
direct HTTP(S) file URLs:

```bash
python -m bioseq.download variants "https://example.org/sample.vcf.gz"
python -m bioseq.download bam "https://example.org/sample.bam"
python -m bioseq.download alignment "https://example.org/sample.cram"
```

VCFs go in `data/variants/`; BAM and CRAM files go in `data/alignments/`.
These commands download existing files—they do not call variants. Direct
downloads require a URL to the actual file; NCBI Datasets does not provide a
VCF-by-variant-accession download through its documented v2 API. Run
`bioseq.pipeline` on a downloaded VCF, BAM, or CRAM to generate its summary
report.

## Next steps

1. Add variant annotation and filtering using an explicitly selected database.
2. Optionally orchestrate downloading, read alignment, and BAM QC in one run.
3. Add biological interpretation while keeping it separate from raw metrics.
