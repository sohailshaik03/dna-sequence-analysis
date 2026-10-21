from fasta_parser import parse_fasta
from dna_analysis import generate_qc_summary


def analyse_fasta(file_path):
    """
    Parse a FASTA file and generate a DNA QC summary
    for each sequence record.

    Returns:
        A list of dictionaries containing each FASTA
        record and its quality-control summary.
    """

    records = parse_fasta(file_path)

    results = []

    for record in records:
        summary = generate_qc_summary(record["sequence"])

        result = {
            "id": record["id"],
            "description": record["description"],
            "qc_summary": summary
        }

        results.append(result)

    return results


def summarise_batch(results):
    """
    Summarise QC results for a batch of FASTA records.

    Args:
        results: Output returned by analyse_fasta().

    Returns:
        A dictionary containing batch-level QC statistics.
    """

    total_records = len(results)

    pass_count = 0
    fail_count = 0
    total_bases = 0

    for result in results:
        summary = result["qc_summary"]

        if summary["qc_status"] == "PASS":
            pass_count += 1
        else:
            fail_count += 1

        total_bases += summary["length"]

    average_length = (
        round(total_bases / total_records, 2)
        if total_records > 0
        else 0
    )

    return {
        "total_records": total_records,
        "pass_count": pass_count,
        "fail_count": fail_count,
        "total_bases": total_bases,
        "average_length": average_length
    }