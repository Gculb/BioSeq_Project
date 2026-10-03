# BioSeq Starter

A small Python bioinformatics pipeline for reasoning about sequence data.

## Planned pipeline

Download data
    ↓
Validate FASTA / FASTQ
    ↓
Basic sequence statistics / QC
    ↓
BLAST
    ↓
Alignment
    ↓
SAM / BAM
    ↓
samtools statistics
    ↓
Biological interpretation

## Current version

The project includes FASTA/FASTQ validation and analysis, FASTQ quality control,
NCBI FASTA and SRA FASTQ download helpers, and JSON/text pipeline reports.
BLAST, alignment, and samtools integration are not part of the current local
pipeline yet.

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

The repository includes small synthetic FASTA and FASTQ examples under
`data/examples/`; input files in `data/raw/` remain ignored and local.

Reports are written to `results/` by default. Override the destination or the
mean-read Phred threshold used for FASTQ low-quality percentages with:

```bash
python -m bioseq.pipeline data/raw/example.fastq \
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

## Suggested development order

1. FASTA parsing
2. FASTQ parsing
3. Sequence statistics
4. FASTQ quality control
5. Data validation
6. Downloading sequence data
7. BLAST integration
8. Alignment / SAM-BAM handling
9. samtools integration
10. Report generation
