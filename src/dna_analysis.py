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



if __name__ == "__main__":
    sequence = "ATGCGTAC"

    print("Sequence:", sequence)
    print("Valid DNA:", validate_dna(sequence))
    print("Length:", calculate_length(sequence))