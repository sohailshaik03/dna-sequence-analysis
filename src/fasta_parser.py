
def parse_fasta(file_path):
    """
    Parse a FASTA file containing one or more sequence records.

    Returns:
        A list of dictionaries containing:
        - id
        - description
        - sequence
    """

    records = []
    sequence_id = None
    description = ""
    sequence_lines = []

    def save_record():
        if not sequence_lines:
            raise ValueError("Invalid FASTA file: empty sequence.")

        records.append({
            "id": sequence_id,
            "description": description,
            "sequence": "".join(sequence_lines).upper()
        })

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                # Save the previous record before starting another.
                if sequence_id is not None:
                    save_record()

                header_content = line[1:].strip()
                parts = header_content.split(maxsplit=1)

                if not parts:
                    raise ValueError(
                        "Invalid FASTA file: header is missing an identifier."
                    )

                sequence_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                sequence_lines = []

            else:
                if sequence_id is None:
                    raise ValueError(
                        "Invalid FASTA file: sequence found before header."
                    )

                sequence_lines.append(line)

    # Save the final record after reaching the end of the file.
    if sequence_id is not None:
        save_record()

    if not records:
        raise ValueError("Invalid FASTA file: file is empty.")

    return records
