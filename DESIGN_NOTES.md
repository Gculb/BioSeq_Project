# BioSeq Design Notes

Use this file to reason about the project before implementing each stage.

## Biological questions

### FASTA
- What kind of biological data is stored in FASTA?
- What does each sequence represent?
- What can sequence length tell us?
- What can GC content tell us?

### FASTQ
- Why does FASTQ contain both sequence and quality information?
- What does a Phred quality score represent?
- Why should low-quality reads potentially be filtered?

### BLAST
- What biological question does BLAST answer?
- What is percent identity?
- What is query coverage?
- What is an E-value?
- Why is similarity not necessarily proof of function?

### Alignment
- What is the reference sequence?
- What does it mean for a read to map?
- What information is represented by SAM/BAM?

### samtools
- Why would we sort a BAM?
- Why would we index a BAM?
- What does mapping rate tell us?
- What does coverage tell us?

## Software questions

- What should happen when an input file is malformed?
- Where should raw data live?
- Where should processed data live?
- What information should be preserved in metadata?
- Which steps should be reproducible?
- Which external programs are required?
