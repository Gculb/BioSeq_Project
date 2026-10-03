def read_fastq(filename):
    reads = []

    with open(filename, "r") as file:
        identifier = file.readline().strip()

        while identifier:
            sequence = file.readline().strip()
            plus = file.readline().strip()
            quality = file.readline().strip()

            if not identifier.startswith("@"):
                raise ValueError("Invalid FASTQ identifier")

            if plus != "+":
                raise ValueError("Invalid FASTQ format")

            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence and quality lengths do not match for {identifier}"
                )

            reads.append({
                "id": identifier[1:],
                "sequence": sequence,
                "quality": quality
            })

            identifier = file.readline().strip()

    return reads

def mean_read_length(reads):
    if not reads:
        return 0

    return sum(len(read["sequence"]) for read in reads) / len(reads)
