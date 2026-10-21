import pytest

from src.fasta_parser import parse_fasta

def test_parse_single_fasta_record(tmp_path):
    fasta_file = tmp_path / "test.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ATGCATGC\n"
    )

    records = parse_fasta(fasta_file)

    assert records == [
        {
            "id": "seq1",
            "description": "",
            "sequence": "ATGCATGC"
        }
    ]





def test_parse_multiline_fasta_sequence(tmp_path):
    fasta_file = tmp_path / "multiline.fasta"
    fasta_file.write_text(
        ">seq2 Example sequence\n"
        "ATGC\n"
        "GTAA\n"
        "CC\n"
    )

    records = parse_fasta(fasta_file)

    assert records[0]["id"] == "seq2"
    assert records[0]["description"] == "Example sequence"
    assert records[0]["sequence"] == "ATGCGTAACC"









def test_parse_multiple_fasta_records(tmp_path):
    fasta_file = tmp_path / "multiple.fasta"

    fasta_file.write_text(
        ">gene_A First gene\n"
        "ATGC\n"
        "GTAA\n"
        ">gene_B Second gene\n"
        "CCGG\n"
        "TTAA\n"
    )

    records = parse_fasta(fasta_file)

    assert len(records) == 2

    assert records[0] == {
        "id": "gene_A",
        "description": "First gene",
        "sequence": "ATGCGTAA"
    }

    assert records[1] == {
        "id": "gene_B",
        "description": "Second gene",
        "sequence": "CCGGTTAA"
    }



def test_empty_fasta_file_raises_error(tmp_path):
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_text("")

    with pytest.raises(ValueError, match="empty"):
        parse_fasta(fasta_file)




def test_header_without_identifier_raises_error(tmp_path):
    fasta_file = tmp_path / "missing_id.fasta"
    fasta_file.write_text(">\nATGC\n")

    with pytest.raises(ValueError, match="identifier"):
        parse_fasta(fasta_file)


def test_sequence_before_header_raises_error(tmp_path):
    fasta_file = tmp_path / "sequence_before_header.fasta"
    fasta_file.write_text("ATGC\n>gene_A\nGTAA\n")

    with pytest.raises(ValueError, match="before header"):
        parse_fasta(fasta_file)




def test_empty_sequence_record_raises_error(tmp_path):
    fasta_file = tmp_path / "empty_sequence.fasta"
    fasta_file.write_text(
        ">gene_A First gene\n"
        ">gene_B Second gene\n"
        "CCGG\n"
    )

    with pytest.raises(ValueError, match="empty sequence"):
        parse_fasta(fasta_file)
