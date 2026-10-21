# DNA Sequence Analysis Toolkit

A Python-based bioinformatics toolkit for DNA sequence validation, nucleotide composition analysis, sequence quality control (QC), FASTA parsing, batch analysis, and complementary sequence generation.

This project is being developed as a practical computational biology portfolio project, with emphasis on scientific correctness, reproducible analysis, automated testing, and maintainable code.

## Project Overview

DNA sequence data can contain canonical nucleotides, ambiguous IUPAC nucleotide codes, alignment gaps, and unsupported characters. Reliable downstream analysis requires explicit sequence validation and clearly defined quality-control rules.

This toolkit provides a modular implementation for analysing individual DNA sequences and processing multiple sequences from FASTA files.

### Current Features

- DNA sequence validation and character classification
- Canonical nucleotide counting and composition
- IUPAC ambiguity-code detection and counting
- Gap and invalid-character detection
- GC, AT, and N-content calculations
- Canonical, ambiguity, gap, and invalid-character percentages
- DNA complement and reverse-complement generation
- Structured per-sequence QC summaries
- FASTA parsing for single and multiple records
- Support for multiline FASTA sequences
- Batch-level statistics across FASTA records
- Command-line interface for FASTA analysis
- Automated tests for analytical functions, parsing, batch analysis, and CLI behaviour

## Scientific Scope

### Accepted Sequence Symbols

| Category | Symbols | Interpretation |
|---|---|---|
| Canonical bases | `A`, `T`, `G`, `C` | Standard DNA nucleotides |
| IUPAC ambiguity codes | `N`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V` | Nucleotide uncertainty |
| Gap | `-` | Alignment gap |
| Invalid | Unsupported characters | Not recognised by the toolkit |

Ambiguity codes are accepted as valid symbols but may affect the suitability of a sequence for a particular downstream analysis. Gaps are supported for sequence data that may originate from alignments.

The FASTA parser reads records and normalises sequence characters to uppercase. DNA-specific validation is performed by the analysis module.

### Quality-Control Policy

The current QC system distinguishes sequence validity from composition metrics.

- `PASS`: the sequence is non-empty and contains no unsupported characters.
- `FAIL`: the sequence is empty or contains one or more unsupported characters.

Ambiguity codes and gaps do not automatically cause failure. The toolkit reports these characteristics so that downstream workflows can apply appropriate, analysis-specific criteria.

**Important limitation:** `PASS` indicates compliance with the toolkit's current character-validation rules. It does not establish biological correctness, contamination-free sequence data, or suitability for every downstream analysis.

## Reported Metrics

The integrated QC summary includes:

- Sequence length
- QC status and validity
- Canonical base counts
- IUPAC ambiguity counts
- Character classification
- Invalid characters
- Canonical, ambiguity, gap, and invalid percentages
- N content
- GC and AT content
- Complement sequence
- Reverse-complement sequence

GC and AT percentages use canonical nucleotides only:

```text
GC% = (G + C) / (A + T + G + C) × 100

AT% = (A + T) / (A + T + G + C) × 100
```

Ambiguous bases and gaps are excluded from these denominators. The calculations are descriptive metrics, not universal biological quality thresholds.

## FASTA Analysis

The parser supports:

- Single or multiple FASTA records
- Multiline sequences
- Sequence identifiers and optional descriptions
- Blank lines
- Uppercase normalisation
- Error reporting for malformed input, including missing identifiers and empty sequence records

Example input:

```text
>seq1 Valid DNA sequence
ATGCATGC
>seq2 DNA with ambiguous bases
ATGCNN
>seq3 DNA with invalid character
ATGX
```

## Command-Line Usage

Run the CLI from the repository root.

Analyse a FASTA file:

```bash
PYTHONPATH=src /usr/bin/python3 -m fasta_cli data/example.fasta
```

The command reports each sequence's identifier, description, length, GC content, and QC status, followed by batch statistics.

Example batch summary:

```text
--- Batch Summary ---
Total records: 3
Passed: 2
Failed: 1
Total bases: 18
Average length: 6.0 bp
```

The example assumes that `data/example.fasta` exists. Replace the path with the location of your own FASTA file.

The CLI reports user-facing errors for file-access problems and malformed FASTA input.

## Project Structure

```text
dna-sequence-analysis/
├── .vscode/
│   └── settings.json
├── data/
│   └── example.fasta
├── src/
│   ├── dna_analysis.py
│   ├── fasta_parser.py
│   ├── fasta_analysis.py
│   └── fasta_cli.py
├── tests/
│   ├── test_dna_analysis.py
│   ├── test_fasta_parser.py
│   ├── test_fasta_analysis.py
│   └── test_fasta_cli.py
├── .gitattributes
├── .gitignore
└── README.md
```

`data/example.fasta` is a local demonstration file. The `.vscode/settings.json` file configures Python analysis to recognise the `src/` directory.

### Module Responsibilities

- `dna_analysis.py`: sequence validation, composition metrics, complement generation, and QC summaries.
- `fasta_parser.py`: FASTA file parsing and record extraction.
- `fasta_analysis.py`: per-record QC analysis and batch statistics.
- `fasta_cli.py`: terminal interface for running the analysis.
- `tests/`: automated verification of the individual modules and integrated workflow.

## Running the Tests

The project uses `pytest`.

Run the complete test suite from the repository root:

```bash
PYTHONPATH=src /usr/bin/python3 -m pytest -v
```

Run a specific test module:

```bash
PYTHONPATH=src /usr/bin/python3 -m pytest tests/test_fasta_parser.py -v
```

The latest verified test suite contains **106 passing tests**. This count reflects the current development state and may change as new features and tests are added.

## Technologies

- Python 3
- pytest
- Git and GitHub
- Visual Studio Code

## Design Principles

1. **Explicit biological assumptions:** sequence-processing rules are defined rather than hidden in calculations.
2. **Modularity:** parsing, analysis, and command-line execution have separate responsibilities.
3. **Validation before interpretation:** unsupported characters are identified rather than silently removed.
4. **Reproducibility:** deterministic functions and automated tests support repeatable results.
5. **Regression testing:** automated tests help detect unintended changes.
6. **Separation of QC and interpretation:** the toolkit reports sequence characteristics without imposing experiment-specific biological thresholds.

## Current Limitations and Future Development

The current implementation focuses on basic DNA sequence QC and FASTA processing. It does not yet implement:

- Open reading frame (ORF) identification
- Codon analysis and translation
- Motif detection
- GC-content visualisation
- Configurable biological QC thresholds
- Large-scale performance benchmarking
- Continuous integration with GitHub Actions
- Formal packaging and installation as a command-line application

These are potential future extensions rather than current features.

## Author

**Sohail**

MSc Bioinformatics

This project is part of a broader portfolio focused on computational biology, bioinformatics, data analysis, and reproducible scientific programming.