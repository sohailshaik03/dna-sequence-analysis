from src.fasta_analysis import analyse_fasta, summarise_batch

def test_analyse_single_fasta_record(tmp_path):
    fasta_file = tmp_path / "single.fasta"

    fasta_file.write_text(
        ">seq1 Example DNA sequence\n"
        "ATGC\n"
    )

    results = analyse_fasta(str(fasta_file))

    assert len(results) == 1
    assert results[0]["id"] == "seq1"
    assert results[0]["description"] == "Example DNA sequence"
    assert results[0]["qc_summary"]["length"] == 4
    assert results[0]["qc_summary"]["gc_content"] == 50.0
    assert results[0]["qc_summary"]["qc_status"] == "PASS"


def test_analyse_multiple_fasta_records(tmp_path):
    fasta_file = tmp_path / "multiple.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ATGC\n"
        ">seq2\n"
        "AAAA\n"
    )

    results = analyse_fasta(str(fasta_file))

    assert len(results) == 2
    assert results[0]["id"] == "seq1"
    assert results[1]["id"] == "seq2"
    assert results[0]["qc_summary"]["gc_content"] == 50.0
    assert results[1]["qc_summary"]["gc_content"] == 0.0


def test_analyse_fasta_record_with_invalid_dna(tmp_path):
    fasta_file = tmp_path / "invalid.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ATGX\n"
    )

    results = analyse_fasta(str(fasta_file))
    summary = results[0]["qc_summary"]

    assert summary["valid"] is False
    assert summary["invalid_bases"] == ["X"]
    assert summary["qc_status"] == "FAIL"
    assert summary["complement"] is None




def test_summarise_batch_counts_and_lengths(tmp_path):
    fasta_file = tmp_path / "batch.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ATGC\n"
        ">seq2\n"
        "AAAAAA\n"
        ">seq3\n"
        "ATGX\n"
    )

    results = analyse_fasta(str(fasta_file))
    summary = summarise_batch(results)

    assert summary["total_records"] == 3
    assert summary["pass_count"] == 2
    assert summary["fail_count"] == 1
    assert summary["total_bases"] == 14
    assert summary["average_length"] == 4.67


def test_summarise_empty_results():
    summary = summarise_batch([])

    assert summary["total_records"] == 0
    assert summary["pass_count"] == 0
    assert summary["fail_count"] == 0
    assert summary["total_bases"] == 0
    assert summary["average_length"] == 0

def test_summarise_batch_counts_and_lengths(tmp_path):
    fasta_file = tmp_path / "batch.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ATGC\n"
        ">seq2\n"
        "AAAAAA\n"
        ">seq3\n"
        "ATGX\n"
    )

    results = analyse_fasta(str(fasta_file))
    summary = summarise_batch(results)

    assert summary["total_records"] == 3
    assert summary["pass_count"] == 2
    assert summary["fail_count"] == 1
    assert summary["total_bases"] == 14
    assert summary["average_length"] == 4.67


def test_summarise_empty_results():
    summary = summarise_batch([])

    assert summary["total_records"] == 0
    assert summary["pass_count"] == 0
    assert summary["fail_count"] == 0
    assert summary["total_bases"] == 0
    assert summary["average_length"] == 0