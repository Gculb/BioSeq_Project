def read_fastq(filename):
    reads = []

    with open(filename, "r", encoding="utf-8") as file_handle:
        lines = [line.rstrip("\r\n") for line in file_handle]

    if len(lines) % 4:
        raise ValueError("FASTQ file ends with an incomplete four-line record.")

    for index in range(0, len(lines), 4):
        identifier, sequence, plus, quality = lines[index : index + 4]
        record_number = index // 4 + 1

        if not identifier.startswith("@") or len(identifier) == 1:
            raise ValueError(f"Invalid FASTQ identifier on record {record_number}.")

        if not plus.startswith("+"):
            raise ValueError(f"Invalid FASTQ separator on record {record_number}.")

        if len(sequence) != len(quality):
            raise ValueError(
                f"Sequence and quality lengths do not match for {identifier}"
            )

        reads.append(
            {"id": identifier[1:], "sequence": sequence, "quality": quality}
        )

    return reads

def mean_read_length(reads):
    if not reads:
        return 0

    return sum(len(read["sequence"]) for read in reads) / len(reads)
