# DNA Sequence Analysis Toolkit

A Python-based bioinformatics toolkit for DNA sequence validation, nucleotide composition analysis, sequence quality control (QC), and complementary sequence generation.

The project is developed as a practical bioinformatics software project, with emphasis on **scientific correctness, reproducible analysis, automated testing, and maintainable code**.

---

## Project Overview

DNA sequence data can contain canonical nucleotides, ambiguous IUPAC nucleotide codes, alignment gaps, and unsupported characters. Reliable downstream analysis therefore requires explicit sequence validation and quality-control rules.

This toolkit provides a modular Python implementation for analysing DNA sequences and generating a structured QC summary.

The current implementation supports:

* DNA sequence validation
* Canonical nucleotide classification
* IUPAC ambiguity-code classification
* Gap detection
* Invalid-character detection
* Nucleotide counting
* GC and AT content calculation
* N-content calculation
* Sequence composition percentages
* Canonical, ambiguous, gap, and invalid percentages
* Complement generation
* Reverse-complement generation
* Integrated DNA sequence QC reports

---

## Scientific Scope

### Accepted sequence symbols

The toolkit distinguishes four character categories:

| Category              | Symbols                                               | Interpretation                |
| --------------------- | ----------------------------------------------------- | ----------------------------- |
| Canonical bases       | `A`, `T`, `G`, `C`                                    | Standard DNA nucleotides      |
| IUPAC ambiguity codes | `N`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V` | Nucleotide uncertainty        |
| Gap                   | `-`                                                   | Alignment gap                 |
| Invalid               | Any unsupported character                             | Not recognised by the toolkit |

IUPAC ambiguity codes are treated as valid sequence symbols, but their presence should be considered when deciding whether a sequence is suitable for a specific downstream analysis.

Gaps are supported because they can occur in aligned sequence data. Their interpretation depends on whether the input represents an alignment or an unaligned biological sequence.

---

## Quality Control Methodology

The QC system separates **sequence validity** from **sequence composition metrics**.

A sequence receives:

* `PASS` when it is non-empty and contains no unsupported characters.
* `FAIL` when it is empty or contains one or more unsupported characters.

The current QC policy does **not** automatically fail sequences because they contain ambiguity codes or gaps. These characteristics are reported separately so that downstream analyses can apply analysis-specific thresholds.

### Reported QC metrics

The integrated QC summary reports:

* Sequence length
* QC status
* Validity
* Canonical base counts
* IUPAC ambiguity counts
* Character classification
* Invalid characters
* Canonical percentage
* Ambiguity percentage
* Gap percentage
* Invalid percentage
* N content
* GC content
* AT content
* Complement sequence
* Reverse-complement sequence

---

## GC and AT Content

GC and AT content are calculated using **canonical nucleotides only**.

For example:

```text
GC% = (G + C) / (A + T + G + C) × 100

AT% = (A + T) / (A + T + G + C) × 100
```

Ambiguous bases and gaps are therefore excluded from the denominator of these two metrics.

This distinction is important because treating ambiguous or alignment characters as canonical nucleotides could produce misleading nucleotide-composition estimates.

---

## Example

Example input:

```text
ATGCNRY-SW?X
```

The toolkit identifies:

```text
Canonical bases: A, T, G, C
Ambiguous bases: N, R, Y, S, W
Gap: -
Invalid characters: ?, X
```

Because unsupported characters are present, the sequence receives:

```text
QC Status: FAIL
```

The invalid characters are reported explicitly rather than silently removed.

---

## Project Structure

```text
dna-sequence-analysis/
│
├── src/
│   └── dna_analysis.py
│
├── tests/
│   └── test_dna_analysis.py
│
├── .gitignore
├── .gitattributes
└── README.md
```

### Source code

`src/dna_analysis.py` contains the sequence-analysis and QC functions.

### Tests

`tests/test_dna_analysis.py` contains the automated test suite covering individual functions and integrated QC behaviour.

---

## Testing

The project uses **pytest** for automated testing.

Run the complete test suite from the repository root:

```bash
/usr/bin/python3 -m pytest
```

The current test suite contains **92 tests** covering:

* Input validation
* Character classification
* Sequence classification
* Base counting
* Ambiguity counting
* GC/AT/N content
* Percentage calculations
* Base composition
* Complement and reverse complement
* QC status determination
* Integrated QC summaries
* Invalid-base detection
* Sequence-length calculation
* QC report generation

All current tests pass.

---

## Technologies

* **Python 3**
* **pytest**
* **Git**
* **GitHub**
* **Visual Studio Code**

---

## Design Principles

The project follows several principles relevant to scientific software development:

1. **Explicit biological assumptions**
   Sequence-processing rules are defined rather than hidden inside calculations.

2. **Modular functions**
   Individual analytical operations are separated into reusable functions.

3. **Validation before interpretation**
   Unsupported sequence characters are identified before downstream interpretation.

4. **Reproducibility**
   Analysis is implemented as deterministic Python functions with automated tests.

5. **Test-driven verification**
   Expected behaviour is explicitly tested to reduce regression risk.

6. **Separation of QC and interpretation**
   The toolkit reports sequence characteristics without imposing experiment-specific biological thresholds.

---

## Current Status

The core DNA sequence QC functionality is implemented and covered by an automated test suite.

The project is actively being developed toward a more complete bioinformatics sequence-analysis toolkit.

---

## Future Development

Planned extensions include:

* FASTA file parsing
* Batch sequence analysis
* Sequence statistics
* ORF identification
* Codon analysis
* Translation
* Motif detection
* GC-content visualisation
* Sequence-quality visualisation
* Command-line interface
* Configuration and reproducible workflows
* Continuous integration with GitHub Actions
* Expanded documentation and examples

---

## Author

**Sohail**

MSc Bioinformatics

This project is part of a broader portfolio focused on computational biology, bioinformatics, data analysis, and reproducible scientific programming.
