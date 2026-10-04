# BioSeq

[![BioSeq CLI](https://github.com/Gculb/BioSeq_Project/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Gculb/BioSeq_Project/actions/workflows/ci.yml)
[![RNA-seq pipeline](https://github.com/Gculb/BioSeq_Project/actions/workflows/rnaseq.yml/badge.svg?branch=main)](https://github.com/Gculb/BioSeq_Project/actions/workflows/rnaseq.yml)
[![HG002 pipeline](https://github.com/Gculb/BioSeq_Project/actions/workflows/hg002.yml/badge.svg?branch=main)](https://github.com/Gculb/BioSeq_Project/actions/workflows/hg002.yml)

Reproducible genomics workflows with Python, R, Nextflow, and containers.

See [ARCHITECTURE.md](./ARCHITECTURE.md) for component and data-flow diagrams,
and [DESIGN_NOTES.md](./DESIGN_NOTES.md) for answers to the biological and
software questions behind the project.

## The problem and what the project demonstrates

Genomics results depend on a chain of tools, inputs, parameters, and runtime
environments. When those steps live in disconnected scripts, it is difficult
to reproduce a result, see where a sample failed, or tell whether a wrapper
changed the behavior of the underlying bioinformatics tools. BioSeq addresses
that workflow problem with containerized, stage-based pipelines and explicit
quality-control and result artifacts. It is a reproducibility and workflow
engineering project; it does not claim to introduce a new variant caller or
improve the accuracy of established tools.

Two pilots make that goal measurable:

| Pilot | Question | Observed result |
| --- | --- | --- |
| HG002 regional variant calling | Do BioSeq's wrappers preserve the results of the same direct BWA, samtools, and GATK commands? | On one 1 Mb chr20 region, both paths had identical GIAB truth-set scores: precision, recall, and F1 were each 0.9962 (1,582 TP, 6 FP, 6 FN). The wrapper run took 302 s vs 328 s and used 561.6 vs 584.5 MiB peak sampled RSS in this single run. |
| Airway smooth-muscle RNA-seq | Can public paired raw reads be taken through QC and quantification to a donor-adjusted biological comparison? | All eight libraries completed fastp and Salmon processing; Salmon mapping rates were 92.3–93.9%. In the DESeq2 model (`~ cell + dex`), no genes met BH-adjusted `p < 0.05` across four donor cell lines. With only four donors, this is inconclusive—not evidence of no treatment effect. |

These are bounded demonstrations, not general performance or biological
claims. The HG002 runtime and memory differences come from a single run and
need repeated measurements to establish a reliable performance difference.
The RNA-seq analysis is an exploratory reanalysis of one study, not an
independent cohort or a machine-learning model; its principal result is that
this small paired dataset did not yield statistically significant genes at
the selected threshold. See the [HG002 results](#hg002-chr20-pilot-results),
[RNA-seq analysis](#raw-read-rna-seq-differential-expression-analysis), and
[pipeline analysis report](./PIPELINE_ANALYSIS.md) for methods, figures, and
limitations.

## Current version

The project includes FASTA/FASTQ validation and analysis, FASTQ quality control,
NCBI FASTA and SRA FASTQ download helpers, JSON/text pipeline reports, an
HG002 regional variant-calling benchmark, and a Nextflow/RNA-seq analysis.
The local `bioseq.pipeline` processes **one input file at a time**. It accepts
FASTA, FASTQ, VCF/VCF.GZ, BAM, and CRAM files and writes format-specific
quality-control/statistics reports. It does not automatically download data,
run BLAST, align reads, or call variants.

| File type | Can download? | Accepted by `bioseq.pipeline`? | Available processing |
| --- | --- | --- | --- |
| FASTA | Yes | Yes | Validate, sequence counts/lengths, GC content, JSON/text report |
| FASTQ | Yes | Yes | Validate, read-length and Phred quality metrics, GC/ambiguous-base metrics, JSON/text report |
| VCF / VCF.GZ | Yes, from a direct URL | Yes | Variant-record, allele, filter, sample-genotype, and QUAL summaries; no annotation or interpretation |
| BAM | Yes, from a direct URL | Yes | `samtools quickcheck`, `flagstat`, and `stats` summaries |
| CRAM | Yes, from a direct URL | Yes | Same samtools summaries; may need a reference FASTA |

The `bioseq.pipeline` analysis command does not download files automatically.
Download and analysis remain separate steps.

## GitHub Actions CLI checks

The `BioSeq CLI` GitHub Actions workflow runs on pushes, pull requests, or
manual dispatch. It installs the Python dependencies, runs the unit-test suite,
then invokes `bioseq.pipeline` on the example FASTA and FASTQ files and checks
that each run produced JSON and text reports. It does not download external
datasets or run the resource-intensive GIAB benchmark.

The `BioSeq RNA-seq` workflow validates Nextflow on pushes and pull requests.
Its full raw-read analysis is manually dispatched because it downloads several
paired-end public sequencing runs. The completed analysis and Nextflow reports
are saved as a workflow artifact.

The `BioSeq HG002 Nextflow benchmark` workflow can be started manually from
GitHub Actions. It runs the regional pilot on a hosted Linux runner; input
preparation is network- and compute-intensive, so it is not part of every push
or pull-request check. Its results and Nextflow execution reports are saved as
the `hg002-nextflow-results` artifact.

## End-to-end HG002 germline pilot

The separate `bioseq.full_benchmark` command runs a small-region short-read
workflow using BioSeq's BWA/samtools wrappers and GATK HaplotypeCaller, then
compares its calls to a direct-command baseline and the GIAB HG002 truth set.
The compared workflows use the same tool versions and parameters; this tests
wrapper parity and overhead, not a claim that BioSeq changes caller accuracy.
The pipeline measures paired FASTQ QC, alignment and duplicate-marking stages,
samtools mapping/coverage metrics, variant metrics, and per-stage elapsed time,
CPU time, sampled peak process-tree RSS, and output size. Reports include
software versions and SHA-256 hashes for inputs. Peak memory is sampled every
50 ms and may miss short-lived peaks; process-tree RSS sums process RSS and can
double-count shared memory.

The included data-preparation command uses the public GIAB HG002 2x250-bp
GRCh38 BAM to extract pairs overlapping a one-megabase chromosome 20 region,
then obtains the matching GIAB truth VCF/confident regions and chromosome
reference plus the Broad GRCh38 Mills/1000G indel known-sites VCF used for
BQSR. The known-sites VCF is a separate resource from the GIAB evaluation truth
set. It uses indexed remote BAM/VCF access and extracts reads from a 5-kb-padded
region, avoiding a whole-file mate search. The source BAM is listed as 122 GB,
and the whole BAM is not intentionally downloaded, though network range
requests may transfer substantial data. Prepared inputs are ignored by Git and
stored locally under `data/benchmark/`.

The same pilot can be run as a Nextflow DAG. Its steps are separate tasks for
input preparation, reference indexing, paired-read QC, BWA-MEM alignment,
SAM/BAM conversion, name sorting, mate fixing, coordinate sorting, duplicate
marking and indexing, BQSR, HaplotypeCaller, VCF QC, and RTG truth evaluation.
Build the pinned tool image, then run from Linux, macOS, or WSL2 with Java 17,
Nextflow 24.10.5, and Docker available:

```bash
docker build -f containers/Dockerfile -t bioseq:latest .
nextflow run hg002.nf \
  -profile docker \
  --interval chr20:10000000-11000000 \
  --threads 2 \
  --outdir results/HG002_nextflow
```

Nextflow publishes QC, variant calls, and RTG summaries under the selected
output directory. The preparation step is deliberately not cached, so each
run retrieves current public resources. This Nextflow DAG executes the direct
tool stages explicitly; the Python `bioseq.full_benchmark` command remains the
separate wrapper-versus-direct parity benchmark.

To rerun the DAG using an already prepared local input directory instead of
downloading/preparing the public data again, add:

```bash
--prepared-dir data/benchmark/giab_hg002_chr20
```

On Windows, start Docker Desktop with its Linux engine first. Build the pinned
tool environment:

```powershell
docker compose build bioseq
```

Prepare the default `chr20:10000000-11000000` pilot (use a new output directory
if you prepare it again):

```powershell
docker compose run --rm bioseq python -m bioseq.giab_data `
  --interval chr20:10000000-11000000 `
  --threads 2 `
  --output-dir data/benchmark/giab_hg002_chr20
```

Run the BioSeq workflow and baseline, then score both call sets against GIAB:

```powershell
docker compose run --rm bioseq python -m bioseq.full_benchmark `
  --read1 data/benchmark/giab_hg002_chr20/HG002_R1.fastq `
  --read2 data/benchmark/giab_hg002_chr20/HG002_R2.fastq `
  --reference-fasta data/benchmark/giab_hg002_chr20/chr20.fa `
  --known-sites-vcf data/benchmark/giab_hg002_chr20/Mills.hg38.region.vcf.gz `
  --reference-sdf data/benchmark/giab_hg002_chr20/chr20.sdf `
  --truth-vcf data/benchmark/giab_hg002_chr20/HG002.truth.region.vcf.gz `
  --confident-regions data/benchmark/giab_hg002_chr20/HG002_chr20_10000000_11000000.confident.bed `
  --interval chr20:10000000-11000000 `
  --threads 2 `
  --output-dir results/HG002_chr20_comparison
```

The comparison JSON contains baseline and BioSeq run manifests, stage-level
performance deltas, and truth-set precision, sensitivity/recall, F1, and
TP/FP/FN counts. Each measured workflow runs once; repeat with a new output
directory and the same inputs to estimate timing variability. Data preparation
and truth-set downloads are reported separately and are not included in the
workflow runtime.

### HG002 chr20 pilot results

The first completed run used the default one-megabase region
`chr20:10000000-11000000`, two threads, and 141,241 paired-end read pairs
extracted from GIAB HG002. The runs used the same pinned tools and parameters;
the baseline invokes the command-line tools directly, while the BioSeq run
invokes them through this project's wrappers.

| Metric | Direct-tool baseline | BioSeq wrappers | Observed difference |
| --- | ---: | ---: | ---: |
| Workflow wall time | 327.71 s | 302.35 s | -25.36 s (-7.7%) |
| Peak sampled process-tree RSS | 584.5 MiB | 561.6 MiB | -22.9 MiB |
| Total intermediate/output size | 575.5 MiB | 575.5 MiB | effectively unchanged |
| Mean alignment depth | 66.80x | 66.80x | equal |
| Bases with at least 1x coverage | 99.9951% | 99.9951% | equal |
| Variant records emitted | 1,749 | 1,749 | equal |

RTG `vcfeval` scored both call sets against the GIAB truth set within its
confident regions:

| Accuracy metric | Baseline | BioSeq |
| --- | ---: | ---: |
| True positives | 1,582 | 1,582 |
| False positives | 6 | 6 |
| False negatives | 6 | 6 |
| Precision | 0.9962 | 0.9962 |
| Sensitivity (recall) | 0.9962 | 0.9962 |
| F1 score | 0.9962 | 0.9962 |

**How to read this:** precision is the proportion of called variants that
match the truth set; sensitivity/recall is the proportion of truth-set variants
recovered; F1 combines the two. A false positive is a call not supported by
truth, and a false negative is a truth variant the workflow missed. The equal
scores are expected here: both runs use the same aligner, caller, tool versions,
and parameters, so this checks that the wrappers preserve behavior. It is not
evidence that the wrappers improve variant accuracy.

The BioSeq run was about 7.7% faster and used 22.9 MiB less peak sampled RSS in
this single measurement. Treat those differences as preliminary, not as a
proven performance improvement: timing can vary with machine load, caching,
and process scheduling. RSS is sampled every 50 ms, sums process RSS (so shared
memory may be counted more than once), and the reported peak is the largest
stage measurement. Output size counts workflow intermediates as well as final
files. The reported coverage is alignment coverage across the requested
interval, not a statement that every base is callable or belongs to the GIAB
confident regions.

The plots below summarize this single regional comparison; the separate
[pipeline analysis report](./PIPELINE_ANALYSIS.md) explains the results and
their limitations.

![HG002 baseline and BioSeq wrapper precision, sensitivity, and F1 parity](./docs/images/hg002-accuracy-parity.svg)

![HG002 baseline and BioSeq wrapper runtime and sampled peak memory](./docs/images/hg002-runtime-memory.svg)

The combined metrics are saved in
[`results/HG002_chr20_comparison/comparison.json`](results/HG002_chr20_comparison/comparison.json);
per-run stage detail is in each implementation's `workflow.json`, and RTG's
truth summaries are under `truth_evaluation/`. These numbers apply only to
this regional pilot, not a full-genome benchmark, a complete GATK Best
Practices workflow (it uses single-sample regional HaplotypeCaller output
rather than GVCF joint genotyping and does not run VQSR), or clinical
validation. BLAST and downloading remain separate utilities; BLAST is not a
read-alignment or variant-calling stage.

## Raw-read RNA-seq differential-expression analysis

The RNA-seq workflow reanalyzes paired raw FASTQ files from the public airway
smooth-muscle study (GEO
[GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778);
SRA study [SRP033351](https://www.ncbi.nlm.nih.gov/sra?term=SRP033351)).
This is a well-known study, but the workflow starts from the original reads,
not the Bioconductor tutorial count matrix. It retrieves the eight untreated
and dexamethasone samples from ENA, runs fastp, quantifies against the
Ensembl 112 GRCh38 transcriptome with Salmon, imports gene-level estimates
with tximport, and performs paired-donor DESeq2 (`~ cell + dex`) in R.

The sample-to-run and treatment mapping is recorded in
[analysis/airway_samples.tsv](./analysis/airway_samples.tsv). Outputs include
per-sample fastp reports, Salmon quantifications, a combined MultiQC report,
adjusted-p-value-ranked gene-level results, a concise dataset-specific
conclusion, MA/volcano/PCA plots, and R session information. MultiQC combines
fastp read-quality metrics with Salmon library/quantification summaries. The
analysis is exploratory, uses four donors from one study, and is not
independent validation or a clinical claim.

In this full raw-read run, no genes passed the Benjamini-Hochberg adjusted
`p < 0.05` threshold after accounting for donor cell line. With only four
donors, this is inconclusive rather than evidence that dexamethasone has no
biological effect. The plots below are generated from this run; the separate
[analysis report](./PIPELINE_ANALYSIS.md) gives the methods and interpretation.

#### Differential-expression visualizations

These plots are generated from the full eight-sample raw-read run and included
here so the analysis can be inspected without downloading the workflow
artifact. The [MultiQC report](./docs/reports/airway-multiqc.html) summarizes
read QC and Salmon quantification across samples.

**PCA of variance-stabilized expression**

![DESeq2 PCA showing airway samples colored by treatment and shaped by donor cell line](./docs/images/airway-pca.png)

**MA plot**

![DESeq2 MA plot for dexamethasone-treated versus untreated samples](./docs/images/airway-ma.png)

**Volcano plot**

![DESeq2 volcano plot of log2 fold change against adjusted p-value](./docs/images/airway-volcano.png)

### Inspecting the results

The manual full-workflow run publishes an `airway-raw-rnaseq-results` artifact.
Open `multiqc/multiqc_report.html` in a browser to inspect QC and Salmon
summaries across samples. The same artifact includes the key DESeq2 outputs:

| Artifact path | What it shows |
| --- | --- |
| `results/analysis_summary.txt` | Dataset-specific conclusion and count of significant genes. |
| `multiqc/multiqc_general_stats.txt` | Machine-readable fastp and Salmon QC summary with one row per sample. |
| `results/ma_plot.png` | Estimated log fold changes across expression levels. |
| `results/volcano_plot.png` | Effect sizes against adjusted p-values. |
| `results/pca_plot.png` | Sample clustering by treatment and donor cell line. |

Build the analysis container from the repository root:

```bash
docker build -f containers/rnaseq.Dockerfile -t bioseq-rnaseq:latest .
```

### Run locally on Windows

Run Nextflow from a WSL2 Linux shell, not from PowerShell. Start Docker Desktop
with its Linux engine and WSL integration enabled for your Linux distribution.
Install Java 17 and Nextflow 24.10.5 inside WSL2, then open the repository in
that shell. The commands below assume the repository is on the Windows C:
drive; adjust the path if yours differs:

```bash
cd "/mnt/c/Users/Gculb/Desktop/BioSeq Project/BioSeq_Project"
```

The RNA-seq run downloads several gigabytes of public reads and needs
additional space for trimmed FASTQs and the Salmon index. Run it from the
repository root:

```bash
docker build -f containers/rnaseq.Dockerfile -t bioseq-rnaseq:latest .
nextflow run main.nf \
  -profile docker \
  --samplesheet analysis/airway_samples.tsv \
  --threads 2 \
  --outdir results/rnaseq_airway_raw
```

When it finishes, view the QC report and DESeq2 plots from WSL2:

```bash
explorer.exe "$(wslpath -w results/rnaseq_airway_raw/multiqc/multiqc_report.html)"
explorer.exe "$(wslpath -w results/rnaseq_airway_raw/results)"
```

The first command opens the combined fastp/Salmon MultiQC report in your
browser; the second opens the folder containing `analysis_summary.txt` and
the MA, volcano, and PCA plots. These files only appear after the full
workflow completes. If Nextflow stops with an error, inspect its terminal
output before trying again.

The HG002 workflow is a separate, optional variant-calling pilot. To run its
staged Nextflow DAG locally, build its tool image and use the prepared inputs
if they are already available:

```bash
docker build -f containers/Dockerfile -t bioseq:latest .
nextflow run hg002.nf \
  -profile docker \
  --prepared-dir data/benchmark/giab_hg002_chr20 \
  --interval chr20:10000000-11000000 \
  --threads 2 \
  --outdir results/HG002_nextflow
```

If that prepared-input directory is absent, omit `--prepared-dir` to let the
workflow retrieve and prepare public inputs. The HG002 DAG produces variant
and truth-scoring reports, not the RNA-seq MultiQC report or DESeq2 plots.

The `docker` profile is configured in `nextflow.config`; the image contains
fastp, Salmon, R, DESeq2, and tximport; a pinned MultiQC container generates
the combined report. On GitHub, use **Actions → BioSeq RNA-seq → Run workflow**
to run the raw-data analysis on a hosted Linux runner, then download its
`airway-raw-rnaseq-results` artifact. The ordinary CI job only validates
workflow syntax to avoid repeatedly downloading large input files on every
push.

## Scope and scaling

The committed HG002 execution is a one-megabase chr20 pilot using two threads;
the RNA-seq example is eight libraries from four donors. These are reproducible
demonstrations, not claims of whole-genome or multi-cohort throughput.

For whole-genome HG002, use full-genome reads and reference inputs rather than
the region-extracted pilot, provision scratch space for large BAM intermediates,
and scatter HaplotypeCaller over non-overlapping intervals before gathering
calls. For multiple individuals, make sample a first-class channel key, align
each sample independently, size CPU/memory/storage per task for the execution
backend, and use GVCF joint genotyping rather than treating samples as one
callset. A real cloud deployment should add a cloud execution profile and
object-store work/results paths; the current hosted-runner examples demonstrate
execution away from a laptop, not a managed cloud deployment.

## Compare variant calls against a truth set

`bioseq.benchmark` compares a baseline caller VCF and a candidate caller VCF
against the same trusted truth VCF using RTG `vcfeval`. It reports true-positive,
false-positive, and false-negative counts plus precision, sensitivity/recall,
F1, and candidate-minus-baseline deltas. This evaluates call-set accuracy; it
does not measure runtime or show that `bioseq.pipeline` itself improved, since
the pipeline currently summarizes VCFs rather than calling variants.

Install RTG Tools and Java, prepare an RTG SDF from the matching reference
assembly, and use the truth set's confident regions where available:

```bash
rtg format -o data/references/GRCh38.sdf data/references/GRCh38.fa
python -m bioseq.benchmark \
  --truth data/benchmark/HG002.truth.vcf.gz \
  --baseline data/benchmark/baseline.vcf.gz \
  --candidate data/benchmark/candidate.vcf.gz \
  --reference data/references/GRCh38.sdf \
  --regions data/benchmark/HG002.confident_regions.bed \
  --output-dir results/HG002_benchmark
```

Use VCFs for the same sample and genome assembly, and do not compare results
outside the regions where the truth set is considered reliable. The current
benchmark command accepts single-sample VCFs and preserves RTG's per-run
`summary.txt` files alongside `benchmark.json`. Output directories are not
overwritten. RTG evaluates PASS variants by default, so keep filtering choices
consistent between baseline and candidate call sets. Results depend on the RTG
version, inputs, reference, and regions; record these when comparing runs.

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
2. Optionally orchestrate downloading, read alignment, variant calling, and QC.
3. Add biological interpretation while keeping it separate from raw metrics.
