def read_fasta(filename):
    sequences = {}

    current_id = None
    current_sequence = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_sequence)

                current_id = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_sequence)

    return sequences


def gc_content(sequence):
    sequence = sequence.upper()

    if not sequence:
        return 0

    gc = sequence.count("G") + sequence.count("C")

    return (gc / len(sequence)) * 100


def sequence_stats(sequences):
    lengths = [len(sequence) for sequence in sequences.values()]

    if not lengths:
        return {
            "count": 0,
            "min_length": 0,
            "max_length": 0,
            "mean_length": 0
        }

    return {
        "count": len(sequences),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": sum(lengths) / len(lengths)
    }
