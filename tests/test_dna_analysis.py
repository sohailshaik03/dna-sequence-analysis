import pytest

from src.dna_analysis import (
    calculate_at_content,
    calculate_ambiguity_percentage,
    calculate_base_composition,
    calculate_canonical_percentage,
    calculate_complement,
    calculate_gap_percentage,
    calculate_invalid_percentage,
    calculate_n_content,
    calculate_percentage,
    calculate_reverse_complement,
    classify_character,
    classify_sequence,
    count_ambiguous_bases,
    count_bases,
    calculate_length,
    determine_qc_status,
    find_invalid_bases,
    gc_content,
    generate_qc_summary,
    print_qc_report,
    validate_dna,
)


def test_validate_dna_accepts_canonical_bases():
    assert validate_dna("ATGC") is True


def test_validate_dna_accepts_lowercase_sequence():
    assert validate_dna("atgc") is True


def test_validate_dna_accepts_iupac_ambiguity_codes():
    assert validate_dna("NRY") is True


def test_validate_dna_accepts_gap():
    assert validate_dna("AT-GC") is True


def test_validate_dna_rejects_invalid_character():
    assert validate_dna("ATGX") is False


def test_validate_dna_rejects_empty_sequence():
    assert validate_dna("") is False



def test_classify_character_canonical():
    assert classify_character("A") == "canonical"


def test_classify_character_ambiguous():
    assert classify_character("N") == "ambiguous"


def test_classify_character_gap():
    assert classify_character("-") == "gap"


def test_classify_character_invalid():
    assert classify_character("X") == "invalid"



def test_classify_sequence_identifies_canonical_bases():
    result = classify_sequence("ATGC")

    assert len(result["canonical"]) == 4
    assert result["canonical"][0] == {
        "position": 1,
        "character": "A"
    }


def test_classify_sequence_identifies_ambiguous_bases():
    result = classify_sequence("NRY")

    assert len(result["ambiguous"]) == 3
    assert result["ambiguous"][0] == {
        "position": 1,
        "character": "N"
    }


def test_classify_sequence_identifies_gaps():
    result = classify_sequence("AT-GC")

    assert len(result["gap"]) == 1
    assert result["gap"][0] == {
        "position": 3,
        "character": "-"
    }


def test_classify_sequence_identifies_invalid_characters():
    result = classify_sequence("ATGX")

    assert len(result["invalid"]) == 1
    assert result["invalid"][0] == {
        "position": 4,
        "character": "X"
    }

def test_count_bases_counts_canonical_bases():
    result = count_bases("AAATTTGGGCCC")

    assert result == {
        "A": 3,
        "T": 3,
        "G": 3,
        "C": 3,
    }


def test_count_bases_handles_lowercase_sequence():
    result = count_bases("aaatttgggccc")

    assert result == {
        "A": 3,
        "T": 3,
        "G": 3,
        "C": 3,
    }


def test_count_bases_ignores_noncanonical_characters():
    result = count_bases("ATGCNRY-SW")

    assert result == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1,
    }


def test_count_ambiguous_bases_counts_iupac_codes():
    result = count_ambiguous_bases("NNRRYYSSWWKKMMBB")

    assert result == {
        "B": 2,
        "D": 0,
        "H": 0,
        "K": 2,
        "M": 2,
        "N": 2,
        "R": 2,
        "S": 2,
        "V": 0,
        "W": 2,
        "Y": 2,
    }


def test_count_ambiguous_bases_handles_lowercase_sequence():
    result = count_ambiguous_bases("nnrryy")

    assert result == {
        "B": 0,
        "D": 0,
        "H": 0,
        "K": 0,
        "M": 0,
        "N": 2,
        "R": 2,
        "S": 0,
        "V": 0,
        "W": 0,
        "Y": 2,
    }


def test_count_ambiguous_bases_ignores_canonical_bases_and_gaps():
    result = count_ambiguous_bases("ATGCNRY-")

    assert result == {
        "B": 0,
        "D": 0,
        "H": 0,
        "K": 0,
        "M": 0,
        "N": 1,
        "R": 1,
        "S": 0,
        "V": 0,
        "W": 0,
        "Y": 1,
    }



def test_gc_content_calculates_percentage():
    assert gc_content("ATGC") == 50.0


def test_gc_content_excludes_ambiguous_bases():
    assert gc_content("ATGCNN") == 50.0


def test_gc_content_excludes_gaps():
    assert gc_content("AT-GC") == 50.0


def test_gc_content_returns_zero_when_no_canonical_bases():
    assert gc_content("NN--") == 0.0


def test_calculate_at_content_calculates_percentage():
    assert calculate_at_content("ATGC") == 50.0


def test_calculate_at_content_excludes_ambiguous_bases():
    assert calculate_at_content("ATGCNN") == 50.0


def test_calculate_at_content_excludes_gaps():
    assert calculate_at_content("AT-GC") == 50.0


def test_calculate_at_content_returns_zero_when_no_canonical_bases():
    assert calculate_at_content("NN--") == 0.0


def test_calculate_n_content_calculates_percentage():
    assert calculate_n_content("ATGCNN") == 33.33


def test_calculate_n_content_handles_lowercase_n():
    assert calculate_n_content("ATGCnn") == 33.33


def test_calculate_n_content_uses_total_sequence_length():
    assert calculate_n_content("NN--") == 50.0


def test_calculate_n_content_returns_zero_for_empty_sequence():
    assert calculate_n_content("") == 0.0


def test_calculate_percentage_calculates_percentage():
    assert calculate_percentage(25, 100) == 25.0


def test_calculate_percentage_rounds_to_two_decimal_places():
    assert calculate_percentage(1, 3) == 33.33


def test_calculate_percentage_returns_zero_for_zero_total():
    assert calculate_percentage(10, 0) == 0.0


def test_calculate_percentage_handles_zero_count():
    assert calculate_percentage(0, 100) == 0.0


def test_calculate_canonical_percentage_calculates_percentage():
    assert calculate_canonical_percentage("ATGC") == 100.0


def test_calculate_canonical_percentage_excludes_ambiguous_bases():
    assert calculate_canonical_percentage("ATGCNN") == 66.67


def test_calculate_canonical_percentage_excludes_gaps():
    assert calculate_canonical_percentage("AT-GC") == 80.0


def test_calculate_canonical_percentage_returns_zero_for_empty_sequence():
    assert calculate_canonical_percentage("") == 0.0


def test_calculate_gap_percentage_calculates_percentage():
    assert calculate_gap_percentage("AT-GC") == 20.0


def test_calculate_gap_percentage_handles_multiple_gaps():
    assert calculate_gap_percentage("A--GC-") == 50.0


def test_calculate_gap_percentage_returns_zero_when_no_gaps():
    assert calculate_gap_percentage("ATGC") == 0.0


def test_calculate_gap_percentage_returns_zero_for_empty_sequence():
    assert calculate_gap_percentage("") == 0.0


def test_calculate_invalid_percentage_calculates_percentage():
    assert calculate_invalid_percentage("ATGX") == 25.0


def test_calculate_invalid_percentage_handles_multiple_invalid_characters():
    assert calculate_invalid_percentage("ATGXX?") == 50.0


def test_calculate_invalid_percentage_returns_zero_when_no_invalid_characters():
    assert calculate_invalid_percentage("ATGCNRY-SW") == 0.0


def test_calculate_invalid_percentage_returns_zero_for_empty_sequence():
    assert calculate_invalid_percentage("") == 0.0




def test_calculate_ambiguity_percentage_calculates_percentage():
    assert calculate_ambiguity_percentage("ATGCNN") == 33.33


def test_calculate_ambiguity_percentage_counts_all_iupac_codes():
    assert calculate_ambiguity_percentage("NRY") == 100.0


def test_calculate_ambiguity_percentage_ignores_canonical_and_gaps():
    assert calculate_ambiguity_percentage("AT-N") == 25.0


def test_calculate_ambiguity_percentage_returns_zero_for_empty_sequence():
    assert calculate_ambiguity_percentage("") == 0.0




def test_calculate_base_composition_returns_canonical_percentages():
    composition = calculate_base_composition("AATTGGCC")

    assert composition["A"] == 25.0
    assert composition["T"] == 25.0
    assert composition["G"] == 25.0
    assert composition["C"] == 25.0


def test_calculate_base_composition_returns_ambiguity_percentages():
    composition = calculate_base_composition("NRY")

    assert composition["N"] == 33.33
    assert composition["R"] == 33.33
    assert composition["Y"] == 33.33


def test_calculate_base_composition_returns_gap_percentage():
    composition = calculate_base_composition("AT-GC")

    assert composition["-"] == 20.0


def test_calculate_base_composition_excludes_invalid_characters():
    composition = calculate_base_composition("ATGX")

    assert "X" not in composition
    assert composition["A"] == 25.0
    assert composition["T"] == 25.0
    assert composition["G"] == 25.0


def test_calculate_base_composition_handles_lowercase():
    composition = calculate_base_composition("aatgcc")

    assert composition["A"] == 33.33
    assert composition["T"] == 16.67
    assert composition["G"] == 16.67
    assert composition["C"] == 33.33


def test_calculate_base_composition_returns_empty_dictionary_for_empty_sequence():
    assert calculate_base_composition("") == {}


def test_calculate_complement_handles_canonical_bases():
    assert calculate_complement("ATGC") == "TACG"


def test_calculate_complement_handles_all_iupac_codes():
    assert calculate_complement("NRYKMBVDHSW") == "NYRMKVBHDSW"


def test_calculate_complement_preserves_gaps():
    assert calculate_complement("AT-GC") == "TA-CG"


def test_calculate_complement_handles_lowercase():
    assert calculate_complement("atgc") == "tacg".upper()


def test_calculate_complement_rejects_invalid_characters():
    with pytest.raises(ValueError, match="Unsupported character"):
        calculate_complement("ATGX")


def test_calculate_complement_rejects_multiple_invalid_characters():
    with pytest.raises(ValueError, match="Unsupported character"):
        calculate_complement("ATG?X")


def test_calculate_reverse_complement_handles_canonical_bases():
    assert calculate_reverse_complement("ATGC") == "GCAT"


def test_calculate_reverse_complement_handles_iupac_codes():
    assert calculate_reverse_complement("NRYKMBVDHSW") == "WSDHBVKMRYN"


def test_calculate_reverse_complement_preserves_gap_position_after_reversal():
    assert calculate_reverse_complement("AT-GC") == "GC-AT"


def test_calculate_reverse_complement_handles_lowercase():
    assert calculate_reverse_complement("atgc") == "GCAT"


def test_calculate_reverse_complement_rejects_invalid_characters():
    with pytest.raises(ValueError, match="Unsupported character"):
        calculate_reverse_complement("ATGX")


def test_calculate_reverse_complement_returns_empty_for_empty_sequence():
    assert calculate_reverse_complement("") == ""


def test_determine_qc_status_passes_valid_sequence():
    summary = {
        "length": 4,
        "invalid_bases": []
    }

    assert determine_qc_status(summary) == "PASS"


def test_determine_qc_status_fails_sequence_with_invalid_bases():
    summary = {
        "length": 5,
        "invalid_bases": [
            {"position": 5, "character": "X"}
        ]
    }

    assert determine_qc_status(summary) == "FAIL"


def test_determine_qc_status_fails_empty_sequence():
    summary = {
        "length": 0,
        "invalid_bases": []
    }

    assert determine_qc_status(summary) == "FAIL"


def test_determine_qc_status_passes_sequence_with_ambiguity_and_gaps():
    summary = {
        "length": 6,
        "invalid_bases": []
    }

    assert determine_qc_status(summary) == "PASS"




def test_generate_qc_summary_handles_valid_sequence():
    summary = generate_qc_summary("ATGC")

    assert summary["sequence"] == "ATGC"
    assert summary["length"] == 4
    assert summary["valid"] is True
    assert summary["qc_status"] == "PASS"


def test_generate_qc_summary_calculates_character_counts():
    summary = generate_qc_summary("ATGCNRY-SW")

    assert summary["character_counts"] == {
        "canonical": 4,
        "ambiguous": 5,
        "gaps": 1,
        "invalid": 0
    }


def test_generate_qc_summary_calculates_base_counts():
    summary = generate_qc_summary("AATTGGCCNRY")

    assert summary["base_counts"] == {
        "A": 2,
        "T": 2,
        "G": 2,
        "C": 2
    }


def test_generate_qc_summary_calculates_percentages():
    summary = generate_qc_summary("ATGCNRY-SW")

    assert summary["percentages"] == {
        "canonical": 40.0,
        "ambiguous": 50.0,
        "gaps": 10.0,
        "invalid": 0.0
    }


def test_generate_qc_summary_calculates_sequence_content():
    summary = generate_qc_summary("ATGCNN")

    assert summary["n_content"] == 33.33
    assert summary["gc_content"] == 50.0
    assert summary["at_content"] == 50.0


def test_generate_qc_summary_handles_invalid_sequence():
    summary = generate_qc_summary("ATGCX?")

    assert summary["valid"] is False
    assert summary["qc_status"] == "FAIL"
    assert summary["invalid_bases"] == ["?", "X"]
    assert summary["complement"] is None
    assert summary["reverse_complement"] is None



def test_generate_qc_summary_includes_classification():
    summary = generate_qc_summary("AT-NX")

    assert summary["classification"] == {
        "canonical": [
            {"position": 1, "character": "A"},
            {"position": 2, "character": "T"},
        ],
        "ambiguous": [
            {"position": 4, "character": "N"},
        ],
        "gap": [
            {"position": 3, "character": "-"},
        ],
        "invalid": [
            {"position": 5, "character": "X"},
        ],
    }




def test_generate_qc_summary_calculates_ambiguous_counts():
    summary = generate_qc_summary("NRYKMBVDHSW")

    assert summary["ambiguous_counts"] == {
        "B": 1,
        "D": 1,
        "H": 1,
        "K": 1,
        "M": 1,
        "N": 1,
        "R": 1,
        "S": 1,
        "V": 1,
        "W": 1,
        "Y": 1,
    }




def test_generate_qc_summary_handles_empty_sequence():
    summary = generate_qc_summary("")

    assert summary["sequence"] == ""
    assert summary["length"] == 0
    assert summary["valid"] is False
    assert summary["qc_status"] == "FAIL"
    assert summary["invalid_bases"] == []
    assert summary["n_content"] == 0.0
    assert summary["gc_content"] == 0.0
    assert summary["at_content"] == 0.0
    assert summary["complement"] == ""
    assert summary["reverse_complement"] == ""





def test_find_invalid_bases_returns_invalid_characters():
    assert find_invalid_bases("ATGCX?") == ["?", "X"]


def test_find_invalid_bases_returns_unique_invalid_characters():
    assert find_invalid_bases("ATXXX??G") == ["?", "X"]


def test_find_invalid_bases_accepts_lowercase_input():
    assert find_invalid_bases("atgcx?") == ["?", "X"]


def test_find_invalid_bases_returns_empty_list_for_valid_sequence():
    assert find_invalid_bases("ATGCNRY-SW") == []



def test_calculate_length_returns_sequence_length():
    assert calculate_length("ATGC") == 4


def test_calculate_length_handles_lowercase_sequence():
    assert calculate_length("atgc") == 4


def test_calculate_length_includes_gaps_and_ambiguous_bases():
    assert calculate_length("AT-NRY") == 6


def test_calculate_length_handles_empty_sequence():
    assert calculate_length("") == 0




def test_print_qc_report_outputs_summary(capsys):
    summary = generate_qc_summary("ATGC")

    print_qc_report(summary)

    captured = capsys.readouterr()

    assert "DNA SEQUENCE QC REPORT" in captured.out
    assert "Sequence:" in captured.out
    assert "ATGC" in captured.out
    assert "Length: 4 bp" in captured.out
    assert "QC Status: PASS" in captured.out
    assert "Valid: True" in captured.out
    assert "Canonical: 4" in captured.out
    assert "Ambiguous: 0" in captured.out
    assert "Gaps: 0" in captured.out
    assert "Invalid: 0" in captured.out
    assert "GC Content: 50.0 %" in captured.out
    assert "AT Content: 50.0 %" in captured.out
    assert "N Content: 0.0 %" in captured.out
    assert "Invalid Characters:" in captured.out
    assert "None" in captured.out
    assert "Complement:" in captured.out
    assert "TACG" in captured.out
    assert "Reverse Complement:" in captured.out
    assert "GCAT" in captured.out