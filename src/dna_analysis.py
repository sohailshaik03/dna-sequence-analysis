def validate_dna(sequence):
    """
    Validate whether a sequence contains only accepted DNA bases.

    Accepted bases:
    A = Adenine
    T = Thymine
    G = Guanine
    C = Cytosine
    N = Unknown/ambiguous nucleotide
    """

    sequence = sequence.upper()

    valid_bases = {"A", "T", "G", "C", "N"}

    return all(base in valid_bases for base in sequence)


def calculate_length(sequence):
    """
    Calculate the length of a DNA sequence.
    """
    return len(sequence)


def count_bases(sequence):
    """
    Count the occurrences of each base in a DNA sequence.

    Returns a dictionary with bases as keys and their counts as values.
    """
    sequence = sequence.upper()
    counts={"A": sequence.count("A"),
             "T": sequence.count("T"),
             "G": sequence.count("G"),
             "C": sequence.count("C"),
             "N": sequence.count("N")}
    return counts


def gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence.

    GC content is the percentage of bases that are either G or C.
    """
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    total_count = len(sequence)

    if total_count == 0:
        return 0.0

    return (gc_count / total_count) * 100

if __name__ == "__main__":
    sequence = "ATGCGTAC"

    print("Sequence:", sequence)
    print("Valid DNA:", validate_dna(sequence))
    print("Length:", calculate_length(sequence))
    print("Base Counts:", count_bases(sequence))
    print("GC Content:", gc_content(sequence))