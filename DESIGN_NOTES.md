# SeqCheckFlow Design Notes

This guide answers the biological and software questions that inform the
current implementation. For component boundaries and data flow, see
[ARCHITECTURE.md](./ARCHITECTURE.md).

## Biological questions

### FASTA

**What kind of biological data is stored in FASTA?**
FASTA stores biological sequences, commonly DNA, RNA, or protein sequences.
Each record starts with a `>` header followed by sequence text.

**What does each sequence represent?**
A record can represent anything from a short read or gene to a transcript,
contig, chromosome, or protein. Its header supplies an identifier and may
include descriptive metadata; the sequence itself does not tell us which
interpretation is correct.

**What can sequence length tell us?**
Length can help check whether a record is plausible for its intended purpose,
compare records, or summarize an assembly. It does not by itself identify a
sequence or establish biological function.

**What can GC content tell us?**
GC content is the fraction of sequence bases that are G or C. It can describe
base composition and help flag unusual regions or technical biases, but it is
not proof of an organism, gene, or function. Interpret it with the sequence
type, organism, and region in mind.

### FASTQ

**Why does FASTQ contain both sequence and quality information?**
The sequence records the called bases; the quality string records the
sequencer's confidence in each corresponding base. Keeping both lets analysis
account for likely sequencing errors rather than treating every base call as
equally certain.

**What does a Phred quality score represent?**
Phred score `Q` is `-10 log10(p)`, where `p` is the estimated probability that
the base call is wrong. For example, Q20 corresponds to an estimated 1% error
probability and Q30 to 0.1%. This is an error-probability scale, not a guarantee
that an individual base is correct.

**Why should low-quality reads potentially be filtered?**
Low-quality reads or bases can introduce false alignments and variant calls.
Filtering or trimming can reduce those errors, but removes data and may bias
coverage, so thresholds should be chosen and recorded for the experiment.
SeqCheckFlow's FASTQ QC currently reports read-quality metrics and the count below a
mean-quality threshold; it does **not** filter or trim reads.

### BLAST

**What biological question does BLAST answer?**
BLAST searches a selected sequence database for local sequence regions similar
to a query. It can identify candidate homologs or related sequences for
follow-up; it does not directly determine biological function.

**What is percent identity?**
The percentage of aligned positions in a selected alignment segment that have
the same residue or nucleotide. It depends on the alignment and is not the
same as similarity, which can also account for conservative amino-acid
substitutions.

**What is query coverage?**
The fraction of the query sequence spanned by the reported alignment. A
high-identity hit covering a small part of the query may be less informative
than a substantial alignment, so identity and coverage should be considered
together.

**What is an E-value?**
The expected number of matches with a score at least this good that could
occur by chance in a database search of the given size. Lower values are
stronger evidence against a random match under the search model; an E-value is
not the probability that a biological hypothesis is true.

**Why is similarity not necessarily proof of function?**
Similar sequences can share ancestry without retaining the same function;
short or low-complexity matches can be misleading; and function depends on
context, structure, expression, and experimental evidence. BLAST results are
candidate evidence that needs independent validation.

### Alignment

**What is the reference sequence?**
A reference is a chosen coordinate-bearing sequence assembly or sequence set
against which reads are compared. Assembly/version and contig naming matter:
the reference must be appropriate for the sample and downstream resources.

**What does it mean for a read to map?**
An aligner found a plausible placement of the read on the reference, allowing
for mismatches and sometimes gaps. Mapping is a computational match, not proof
that the placement is biologically correct; mapping quality estimates
placement confidence.

**What information is represented by SAM/BAM?**
SAM is a text alignment format and BAM is its compressed binary counterpart.
They can record read names, flags, reference and position, mapping quality,
CIGAR alignment operations, mate information, sequence, base qualities, and
optional tags such as read groups. BAM is generally more compact and faster
for programmatic access.

### samtools

**Why would we sort a BAM?**
Coordinate sorting orders alignments by reference and position, which is
required by many tools and makes index-based regional access possible.
Name-sorting instead groups alignments by read name, useful for operations
such as pairing reads or running `samtools fixmate`. The correct order depends
on the next tool.

**Why would we index a BAM?**
An index maps genomic regions to file blocks so tools can retrieve alignments
from a region without scanning the entire file. Coordinate sorting is normally
required first.

**What does mapping rate tell us?**
It is the fraction of reads that the aligner reports as mapped under the
chosen filters and reference. It can help detect poor read quality, an
inappropriate reference, contamination, or sample/library issues, but a high
rate does not guarantee correct placements or variant calls.

**What does coverage tell us?**
Depth is the number of reads covering a position; breadth is the fraction of
positions covered to a chosen depth. Average depth can conceal uncovered
regions, unevenness, duplicates, and low-quality alignments. Coverage supports
assessment of callability, but does not alone establish that a variant call is
correct. The benchmark reports mean depth and the fraction of requested-region
bases with at least one aligned base.

## Software questions

**What should happen when an input file is malformed?**
SeqCheckFlow should reject it, identify the format and validation problem, and stop
before presenting success-shaped analysis output. The CLI reports an error
and exits unsuccessfully. FASTA/FASTQ/VCF validators check their supported
structural rules; BAM/CRAM integrity checks use `samtools quickcheck`. Passing
format validation does not guarantee biological correctness.

**Where should raw data live?**
Use `data/raw/` for source inputs that should remain unchanged. Raw data is
ignored by Git; do not commit sensitive or large sample data. Download helpers
can also save to `data/` by default, with VCFs under `data/variants/` and
BAM/CRAM under `data/alignments/`.

**Where should processed data live?**
Keep derived data under `data/processed/` when it needs to be retained as an
intermediate. QC and benchmark reports belong in `results/`. The regional GIAB
benchmark inputs live under `data/benchmark/`. Generated/large data directories
are ignored by Git; keep the small directory placeholders as needed.

**What information should be preserved in metadata?**
At minimum preserve sample identity (without exposing private identifiers),
reference assembly and contig/region, input paths and checksums, tool names and
versions, command parameters, filters, known-sites and truth-set versions,
workflow implementation, and run date/environment. For benchmark reproducibility,
also preserve per-stage wall time, CPU/RSS measurements, output sizes, and
truth-evaluation definitions. SeqCheckFlow's full benchmark records many of these
fields and hashes the reads, reference, known-sites file, truth VCF, confident
regions, and reference SDF.

**Which steps should be reproducible?**
All transformations from source inputs to reported metrics should be
re-runnable: data selection and region definition, reference/resource
preparation, validation, alignment, sorting/duplicate handling, BQSR, variant
calling, and truth-set scoring. Pin tool versions and parameters, keep inputs
immutable, record checksums and environment details, and use the same
confident-region/filter definitions when comparing runs. Identical settings
improve repeatability but do not guarantee bit-for-bit-identical timing or
outputs across machines.

**Which external programs are required?**
Requirements depend on the selected path:

- Core Python analysis uses Python plus dependencies in `requirements.txt`
  (currently Biopython and psutil).
- BAM/CRAM validation and summaries require `samtools`.
- Read alignment requires `bwa` (or `minimap2` for the standalone alignment
  wrapper); the integrated benchmark uses BWA.
- The end-to-end germline benchmark additionally uses GATK, `bcftools`, and
  RTG Tools. The supplied Docker environment pins these versions.
- SRA-to-FASTQ download requires NCBI SRA Toolkit's `fasterq-dump`.
- BLAST is submitted to NCBI's remote service through Biopython; it requires
  network access rather than a local BLAST executable.

## Current scope reminder

`seqcheckflow.pipeline` validates and summarizes one supplied FASTA, FASTQ, VCF,
BAM, or CRAM at a time. It does not automatically download, align, or call
variants. The separate HG002 benchmark connects paired FASTQ alignment to
variant calling and truth evaluation for a defined small region. See
[ARCHITECTURE.md](./ARCHITECTURE.md) for the component map and these workflow
boundaries.
