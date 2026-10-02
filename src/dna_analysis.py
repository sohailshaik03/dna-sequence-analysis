

CANONICAL_BASES = {"A", "T", "G", "C"}

IUPAC_AMBIGUOUS_BASES = {
    "N", "R", "Y", "S", "W",
    "K", "M", "B", "D", "H", "V"
}

GAP_SYMBOLS = {"-"}


COMPLEMENT_MAP = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",

    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",

    "-": "-"
}


def classify_character(character):
    """
    Classify a single DNA sequence character.

    Returns:
        "canonical"  -> A, T, G, C
        "ambiguous"  -> IUPAC ambiguity code
        "gap"        -> alignment gap
        "invalid"    -> unsupported character
    """

    if character in CANONICAL_BASES:
        return "canonical"

    if character in IUPAC_AMBIGUOUS_BASES:
        return "ambiguous"

    if character in GAP_SYMBOLS:
        return "gap"

    return "invalid"





def classify_sequence(sequence):
    """
    Classify every character in a DNA sequence.

    Returns a dictionary containing the positions and characters
    belonging to each category.

    Positions are 1-based to make the output easier to interpret
    in biological reports.
    """

    sequence = sequence.upper()

    classification = {
        "canonical": [],
        "ambiguous": [],
        "gap": [],
        "invalid": []
    }

    for position, character in enumerate(sequence, start=1):

        category = classify_character(character)

        classification[category].append({
            "position": position,
            "character": character
        })

    return classification







def validate_dna(sequence):
    """
    Validate whether a sequence contains only accepted
    DNA-related symbols.
    """

    sequence = sequence.upper()

    if sequence == "":
        return False

    return all(
        classify_character(character) != "invalid"
        for character in sequence
    )





def find_invalid_bases(sequence):
    """
    Identify unsupported characters in a DNA sequence.
    """

    sequence = sequence.upper()

    return sorted(
        set(
            character
            for character in sequence
            if classify_character(character) == "invalid"
        )
    )




def calculate_length(sequence):
    """
    Calculate the length of a DNA sequence.
    """
    return len(sequence)





def count_bases(sequence):
    """
    Count canonical DNA bases.

    Returns counts for:
        A, T, G and C.
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    counts = {
        "A": 0,
        "T": 0,
        "G": 0,
        "C": 0
    }

    for item in classification["canonical"]:
        base = item["character"]
        counts[base] += 1

    return counts





def count_ambiguous_bases(sequence):
    """
    Count IUPAC ambiguity codes in a DNA sequence.

    Returns counts for:
        N, R, Y, S, W, K, M, B, D, H, V
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    counts = {
        base: 0
        for base in sorted(IUPAC_AMBIGUOUS_BASES)
    }

    for item in classification["ambiguous"]:
        base = item["character"]
        counts[base] += 1

    return counts




def gc_content(sequence):
    """
    Calculate GC content among canonical DNA bases only.

    Formula:

        GC% = (G + C) / (A + T + G + C) * 100

    Ambiguous bases, gaps and invalid characters are excluded
    from the denominator.
    """

    counts = count_bases(sequence)

    canonical_count = sum(counts.values())

    gc_count = counts["G"] + counts["C"]

    if canonical_count == 0:
        return 0.0

    return round((gc_count / canonical_count) * 100, 2)













def calculate_n_content(sequence):
    """
    Calculate the percentage of sequence positions
    represented specifically by N.

    N represents an unknown nucleotide.
    """

    sequence = sequence.upper()

    n_count = sequence.count("N")
    total_length = len(sequence)

    if total_length == 0:
        return 0.0

    return round(
        (n_count / total_length) * 100,
        2
    )


def calculate_at_content(sequence):
    """
    Calculate AT content among canonical DNA bases only.

    Formula:

        AT% = (A + T) / (A + T + G + C) * 100
    """

    counts = count_bases(sequence)

    canonical_count = sum(counts.values())

    at_count = counts["A"] + counts["T"]

    if canonical_count == 0:
        return 0.0

    return round((at_count / canonical_count) * 100, 2)



def calculate_percentage(count, total):
    """
    Calculate a percentage safely.

    Returns 0.0 when the denominator is zero.
    """

    if total == 0:
        return 0.0

    return round((count / total) * 100, 2)





def calculate_canonical_percentage(sequence):
    """
    Calculate the percentage of sequence positions
    represented by canonical bases A, T, G and C.
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    canonical_count = len(classification["canonical"])
    total_length = len(sequence)

    return calculate_percentage(
        canonical_count,
        total_length
    )




def calculate_gap_percentage(sequence):
    """
    Calculate the percentage of sequence positions
    represented by alignment gaps.
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    gap_count = len(classification["gap"])
    total_length = len(sequence)

    return calculate_percentage(
        gap_count,
        total_length
    )



def calculate_invalid_percentage(sequence):
    """
    Calculate the percentage of sequence positions
    containing unsupported characters.
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    invalid_count = len(classification["invalid"])
    total_length = len(sequence)

    return calculate_percentage(
        invalid_count,
        total_length
    )



def calculate_ambiguity_percentage(sequence):
    """
    Calculate the percentage of sequence positions represented
    by IUPAC ambiguity codes.

    Ambiguous bases include:
        N, R, Y, S, W, K, M, B, D, H, V

    Formula:

        Ambiguity% =
        ambiguous bases / total sequence length * 100
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    ambiguous_count = len(classification["ambiguous"])
    total_length = len(sequence)

    if total_length == 0:
        return 0.0

    return round(
        (ambiguous_count / total_length) * 100,
        2
    )





def calculate_base_composition(sequence):
    """
    Calculate the percentage composition of every recognised
    DNA character category.

    Percentages are calculated relative to the total sequence length.

    Categories include:
        Canonical bases: A, T, G, C
        IUPAC ambiguity codes: N, R, Y, S, W, K, M, B, D, H, V
        Gaps: -
        Invalid characters are excluded from the composition dictionary.
    """

    sequence = sequence.upper()

    classification = classify_sequence(sequence)

    total_length = len(sequence)

    if total_length == 0:
        return {}

    composition = {}

    # Canonical bases
    for base in sorted(CANONICAL_BASES):

        count = sum(
            1
            for item in classification["canonical"]
            if item["character"] == base
        )

        composition[base] = calculate_percentage(
            count,
            total_length
        )

    # IUPAC ambiguity codes
    for base in sorted(IUPAC_AMBIGUOUS_BASES):

        count = sum(
            1
            for item in classification["ambiguous"]
            if item["character"] == base
        )

        composition[base] = calculate_percentage(
            count,
            total_length
        )

    # Alignment gaps
    gap_count = len(classification["gap"])

    composition["-"] = calculate_percentage(
        gap_count,
        total_length
    )

    return composition











def generate_qc_summary(sequence):
    """
    Generate a complete DNA sequence quality-control summary.

    The sequence is classified once. All QC counts and percentages
    are derived from this single classification.

    This avoids repeatedly scanning and classifying the same sequence.
    """

    sequence = sequence.upper()

    # Classify the sequence once
    classification = classify_sequence(sequence)

    # Count each category
    canonical_count = len(classification["canonical"])
    ambiguous_count = len(classification["ambiguous"])
    gap_count = len(classification["gap"])
    invalid_count = len(classification["invalid"])

    # Total sequence length
    total_length = len(sequence)

    # Identify unique invalid characters
    invalid_bases = sorted(
        set(
            item["character"]
            for item in classification["invalid"]
        )
    )

    # Count individual canonical bases
    base_counts = {
        "A": 0,
        "T": 0,
        "G": 0,
        "C": 0
    }

    for item in classification["canonical"]:
        base = item["character"]
        base_counts[base] += 1

    # Count individual IUPAC ambiguity codes
    ambiguous_counts = {
        base: 0
        for base in sorted(IUPAC_AMBIGUOUS_BASES)
    }

    for item in classification["ambiguous"]:
        base = item["character"]
        ambiguous_counts[base] += 1

    # Calculate percentages
    canonical_percentage = calculate_percentage(
        canonical_count,
        total_length
    )

    ambiguous_percentage = calculate_percentage(
        ambiguous_count,
        total_length
    )

    gap_percentage = calculate_percentage(
        gap_count,
        total_length
    )

    invalid_percentage = calculate_percentage(
        invalid_count,
        total_length
    )

    n_count = ambiguous_counts["N"]

    n_percentage = calculate_percentage(
        n_count,
        total_length
    )

    # GC and AT percentages are calculated only
    # among canonical A/T/G/C bases.
    canonical_total = sum(base_counts.values())

    gc_percentage = calculate_percentage(
        base_counts["G"] + base_counts["C"],
        canonical_total
    )

    at_percentage = calculate_percentage(
        base_counts["A"] + base_counts["T"],
        canonical_total
    )


    # Calculate complement information only when
    # the sequence contains no unsupported characters.
    if invalid_count == 0:

        complement = calculate_complement(sequence)

        reverse_complement = calculate_reverse_complement(sequence)

    else:

        complement = None

        reverse_complement = None

    summary = {
        "sequence": sequence,

        "length": total_length,

        "valid": (
            total_length > 0
            and invalid_count == 0
        ),

        "invalid_bases": invalid_bases,

        "base_counts": base_counts,

        "ambiguous_counts": ambiguous_counts,

        "character_counts": {
            "canonical": canonical_count,
            "ambiguous": ambiguous_count,
            "gaps": gap_count,
            "invalid": invalid_count
        },

        "percentages": {
            "canonical": canonical_percentage,
            "ambiguous": ambiguous_percentage,
            "gaps": gap_percentage,
            "invalid": invalid_percentage
        },

        "n_content": n_percentage,

        "gc_content": gc_percentage,

        "at_content": at_percentage,

        "complement": complement,

        "reverse_complement": reverse_complement,

        "classification": classification
    }
    summary["qc_status"] = determine_qc_status(summary)
    return summary










def determine_qc_status(summary):
    """
    Determine the validation status of a DNA sequence.

    This function evaluates whether the sequence contains
    unsupported characters.

    It does not impose biological quality thresholds such as
    maximum N content or minimum GC content. Those thresholds
    depend on the specific experiment and analysis pipeline.
    """

    if summary["length"] == 0:
        return "FAIL"

    if summary["invalid_bases"]:
        return "FAIL"

    return "PASS"




def print_qc_report(summary):
    """
    Print a human-readable DNA sequence quality-control report.
    """

    print("\n========== DNA SEQUENCE QC REPORT ==========")

    print("\nSequence:")
    print(" ", summary["sequence"])

    print("\nSequence Statistics:")
    print("  Length:", summary["length"], "bp")
    print("  QC Status:", summary["qc_status"])
    print("  Valid:", summary["valid"])

    print("\nCharacter Counts:")
    print("  Canonical:", summary["character_counts"]["canonical"])
    print("  Ambiguous:", summary["character_counts"]["ambiguous"])
    print("  Gaps:", summary["character_counts"]["gaps"])
    print("  Invalid:", summary["character_counts"]["invalid"])

    print("\nPercentages:")
    print("  Canonical:", summary["percentages"]["canonical"], "%")
    print("  Ambiguous:", summary["percentages"]["ambiguous"], "%")
    print("  Gaps:", summary["percentages"]["gaps"], "%")
    print("  Invalid:", summary["percentages"]["invalid"], "%")

    print("\nCanonical Base Counts:")
    for base, count in summary["base_counts"].items():
        print(f"  {base}: {count}")

    print("\nAmbiguous Base Counts:")
    for base, count in summary["ambiguous_counts"].items():
        print(f"  {base}: {count}")

    print("\nSequence Composition:")
    print("  GC Content:", summary["gc_content"], "%")
    print("  AT Content:", summary["at_content"], "%")
    print("  N Content:", summary["n_content"], "%")

    print("\nInvalid Characters:")

    if summary["invalid_bases"]:
        print(" ", ", ".join(summary["invalid_bases"]))
    else:
        print("  None")

    print("\nComplement:")
    print(" ", summary["complement"])

    print("\nReverse Complement:")
    print(" ", summary["reverse_complement"])

    print("\n============================================")






def calculate_complement(sequence):
    """
    Calculate the DNA complement of a sequence.

    Supports:
        - Canonical DNA bases
        - IUPAC ambiguity codes
        - Alignment gaps

    Raises:
        ValueError: If the sequence contains an unsupported character.
    """

    sequence = sequence.upper()

    invalid_characters = [
        character
        for character in sequence
        if character not in COMPLEMENT_MAP
    ]

    if invalid_characters:
        invalid_characters = sorted(set(invalid_characters))

        raise ValueError(
            f"Unsupported character(s) for complement: "
            f"{', '.join(invalid_characters)}"
        )

    return "".join(
        COMPLEMENT_MAP[character]
        for character in sequence
    )


def calculate_reverse_complement(sequence):
    """
    Calculate the reverse complement of a DNA sequence.

    The reverse complement is obtained by first finding the complement
    and then reversing the resulting sequence.
    """

    complement_sequence = calculate_complement(sequence)
    return complement_sequence[::-1]



if __name__ == "__main__":

    sequence = "ATGCNRY-SW"

    summary = generate_qc_summary(sequence)

    print_qc_report(summary)