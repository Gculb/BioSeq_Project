from collections import Counter

from .fastq import read_fastq


def analyze_fastq_quality(fastq_file, minimum_mean_quality=20):
    """Calculate read-length, base-quality, GC, and ambiguity metrics.

    Quality characters are interpreted as Phred+33. A read is low quality when
    its mean per-base Phred score is below ``minimum_mean_quality``. GC percent
    uses all sequence bases as its denominator; ambiguous bases count bases
    other than A, C, G, and T.
    """
    if minimum_mean_quality < 0:
        raise ValueError("minimum_mean_quality must not be negative.")

    reads = read_fastq(fastq_file)
    read_lengths = Counter()
    quality_distribution = Counter()
    total_bases = 0
    gc_bases = 0
    ambiguous_bases = 0
    low_quality_reads = 0
    total_quality = 0

    for read in reads:
        sequence = read["sequence"].upper()
        quality_scores = [ord(character) - 33 for character in read["quality"]]
        read_lengths[len(sequence)] += 1
        total_bases += len(sequence)
        gc_bases += sequence.count("G") + sequence.count("C")
        ambiguous_bases += sum(base not in "ACGT" for base in sequence)
        quality_distribution.update(quality_scores)
        total_quality += sum(quality_scores)

        mean_read_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
        if mean_read_quality < minimum_mean_quality:
            low_quality_reads += 1

    read_count = len(reads)
    return {
        "read_count": read_count,
        "total_bases": total_bases,
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "mean_read_length": total_bases / read_count if read_count else 0,
        "mean_quality": total_quality / total_bases if total_bases else 0,
        "quality_distribution": dict(sorted(quality_distribution.items())),
        "gc_content": (gc_bases / total_bases * 100) if total_bases else 0,
        "ambiguous_bases": ambiguous_bases,
        "low_quality_reads": low_quality_reads,
        "low_quality_percentage": (low_quality_reads / read_count * 100) if read_count else 0,
        "minimum_mean_quality": minimum_mean_quality,
    }
