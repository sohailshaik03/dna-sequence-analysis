
import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run_cli(fasta_file):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")

    return subprocess.run(
        [
            sys.executable,
            "-m",
            "fasta_cli",
            str(fasta_file),
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )


def test_cli_reports_sequence_and_batch_summary(tmp_path):
    fasta_file = tmp_path / "example.fasta"
    fasta_file.write_text(
        ">seq1 Valid DNA\n"
        "ATGCATGC\n"
        ">seq2 Ambiguous DNA\n"
        "ATGCNN\n"
        ">seq3 Invalid DNA\n"
        "ATGX\n"
    )

    result = run_cli(fasta_file)

    assert result.returncode == 0
    assert "Sequence ID: seq1" in result.stdout
    assert "QC status: PASS" in result.stdout
    assert "Sequence ID: seq3" in result.stdout
    assert "QC status: FAIL" in result.stdout
    assert "Total records: 3" in result.stdout
    assert "Passed: 2" in result.stdout
    assert "Failed: 1" in result.stdout
    assert "Total bases: 18" in result.stdout
    assert "Average length: 6.0 bp" in result.stdout


def test_cli_reports_missing_file(tmp_path):
    missing_file = tmp_path / "missing.fasta"

    result = run_cli(missing_file)

    assert result.returncode == 2
    assert "No such file or directory" in result.stderr
    assert "Traceback" not in result.stderr

def test_cli_reports_malformed_fasta(tmp_path):
    fasta_file = tmp_path / "malformed.fasta"
    fasta_file.write_text("ATGC\n")

    result = run_cli(fasta_file)

    assert result.returncode == 2
    assert "sequence found before header" in result.stderr
    assert "Traceback" not in result.stderr