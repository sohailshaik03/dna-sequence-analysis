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




def find_invalid_bases(sequence):
    """
    Identify characters that are not accepted DNA bases.
    """

    sequence = sequence.upper()

    valid_bases = {"A", "T", "G", "C", "N"}

    invalid_bases = sorted(
        set(base for base in sequence if base not in valid_bases)
    )

    return invalid_bases









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

    return round((gc_count / total_count) * 100, 2)




def calculate_n_content(sequence):
    """
    Calculate the N content of a DNA sequence.

    N content is the percentage of bases that are N (unknown/ambiguous).
    """
    sequence = sequence.upper()
    n_count = sequence.count("N")
    total_count = len(sequence)

    if total_count == 0:
        return 0.0

    return round((n_count / total_count) * 100, 2)


def calculate_at_content(sequence):
    """
    Calculate the AT content of a DNA sequence.

    AT content is the percentage of bases that are either A or T.
    """
    sequence = sequence.upper()
    at_count = sequence.count("A") + sequence.count("T")
    total_count = len(sequence)

    if total_count == 0:
        return 0.0

    return round((at_count / total_count) * 100, 2)




def calculate_base_composition(sequence):
    """
    Calculate the base composition of a DNA sequence.

    Returns a dictionary with bases as keys and their percentage composition as values.
    """
    sequence = sequence.upper()
    total_count = len(sequence)

    if total_count == 0:
        return {"A": 0.0, "T": 0.0, "G": 0.0, "C": 0.0, "N": 0.0}

    composition = {
        "A": round((sequence.count("A") / total_count) * 100, 2),
        "T": round((sequence.count("T") / total_count) * 100, 2),
        "G": round((sequence.count("G") / total_count) * 100, 2),
        "C": round((sequence.count("C") / total_count) * 100, 2),
        "N": round((sequence.count("N") / total_count) * 100, 2)
    }

    return composition



def calculate_known_base_percentage(sequence):
    """
    Calculate the percentage of known DNA bases
    (A, T, G, C) in a sequence.
    """
    sequence = sequence.upper()

    length = len(sequence)

    if length == 0:
        return 0

    known_count = (
        sequence.count("A")
        + sequence.count("T")
        + sequence.count("G")
        + sequence.count("C")
    )

    known_percentage = round((known_count / length) * 100, 2)

    return known_percentage




def generate_qc_summary(sequence):
    """
    Generate a basic quality-control summary
    for a DNA sequence.
    """

    sequence = sequence.upper()

    valid = validate_dna(sequence)
    invalid_bases = find_invalid_bases(sequence)

    summary = {
        "sequence": sequence,
        "valid": valid,
        "invalid_bases": invalid_bases,
        "length": calculate_length(sequence)
    }

    if valid:
        summary["base_counts"] = count_bases(sequence)
        summary["gc_content"] = gc_content(sequence)
        summary["at_content"] = calculate_at_content(sequence)
        summary["n_content"] = calculate_n_content(sequence)
        summary["known_base_percentage"] = calculate_known_base_percentage(sequence)
        summary["base_percentages"] = calculate_base_composition(sequence)

        summary["complement"] = calculate_complement(sequence)
        summary["reverse_complement"] = calculate_reverse_complement(sequence)
    summary["qc_status"] = determine_qc_status(summary)
    return summary







def determine_qc_status(summary):
    """
    Determine the quality-control status of a DNA sequence.

    Returns:
        FAIL if the sequence contains invalid bases.
        PASS WITH WARNING if the sequence contains ambiguous bases (N).
        PASS if the sequence is valid and contains no ambiguous bases.
    """

    if not summary["valid"]:
        return "FAIL"

    if summary["n_content"] > 0:
        return "PASS WITH WARNING"

    return "PASS"





def print_qc_report(summary):
    """
    Print a human-readable quality-control report.
    """

    print("\n========== DNA SEQUENCE QC REPORT ==========")

    print("Sequence:", summary["sequence"])
    print("Valid DNA:", summary["valid"])
    print("Length:", summary["length"], "bp")
    print("QC Status:", summary["qc_status"])
    if not summary["valid"]:
        print("\nValidation:")
        print("  Status: FAIL")
        print("  Invalid bases:", ", ".join(summary["invalid_bases"]))

        print("\nQC Status: FAIL")
        print("============================================")

        return

    print("\nValidation:")
    print("  Status: PASS")

    print("\nBase Counts:")
    for base, count in summary["base_counts"].items():
        print(f"  {base}: {count}")

    print("\nBase Composition:")
    for base, percentage in summary["base_percentages"].items():
        print(f"  {base}: {percentage}%")

    print("\nSequence Statistics:")
    print("  GC Content:", summary["gc_content"], "%")
    print("  AT Content:", summary["at_content"], "%")
    print("  N Content:", summary["n_content"], "%")
    print("  Known Bases:", summary["known_base_percentage"], "%")

    print("\nComplement:", summary["complement"])
    print("Reverse Complement:", summary["reverse_complement"])

    print("\nQC Status:", summary["qc_status"])
    print("============================================")








def calculate_complement(sequence):
    """
    Calculate the complement of a DNA sequence.

    A <-> T
    G <-> C
    N -> N (remains unchanged)
    """

    sequence= sequence.upper()

    complement={
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }

    return "".join(complement[base] for base in sequence)




def calculate_reverse_complement(sequence):
    """
    Calculate the reverse complement of a DNA sequence.

    The reverse complement is obtained by first finding the complement
    and then reversing the resulting sequence.
    """

    complement_sequence = calculate_complement(sequence)
    return complement_sequence[::-1]







    
if __name__ == "__main__":
    sequence = "ATGCNA"

   
    summary = generate_qc_summary(sequence)

    print_qc_report(summary)
    
