import argparse

from .fasta import read_fasta, sequence_stats
from .fastq import read_fastq, mean_read_length


def analyze_fasta(filename):
    sequences = read_fasta(filename)
    stats = sequence_stats(sequences)

    print("\nFASTA Analysis")
    print("================")
    print(f"Sequences:       {stats['count']}")
    print(f"Minimum length:  {stats['min_length']} bp")
    print(f"Maximum length:  {stats['max_length']} bp")
    print(f"Mean length:     {stats['mean_length']:.2f} bp")


def analyze_fastq(filename):
    reads = read_fastq(filename)

    print("\nFASTQ Analysis")
    print("================")
    print(f"Reads:            {len(reads)}")
    print(f"Mean read length: {mean_read_length(reads):.2f} bp")


def main():
    parser = argparse.ArgumentParser(
        description="BioSeq: a small bioinformatics sequence analysis toolkit"
    )

    parser.add_argument("filename", help="FASTA or FASTQ file to analyze")

    args = parser.parse_args()

    if args.filename.endswith((".fasta", ".fa", ".fna")):
        analyze_fasta(args.filename)
    elif args.filename.endswith((".fastq", ".fq")):
        analyze_fastq(args.filename)
    else:
        raise ValueError("Unsupported file type")


if __name__ == "__main__":
    main()
