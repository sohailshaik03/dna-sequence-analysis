
import argparse

from fasta_analysis import analyse_fasta, summarise_batch


def main():
    parser = argparse.ArgumentParser(
        description=(
    "Analyse DNA sequences in a FASTA file, report "
    "per-sequence quality control (QC) results, "
    "and summarise the batch."
)
    )

    parser.add_argument(
        "fasta_file",
        help="Path to the input FASTA file."
    )


    args = parser.parse_args()

    try:
        results = analyse_fasta(args.fasta_file)
    except (OSError, ValueError) as error:
        parser.error(str(error))


    for result in results:
        summary = result["qc_summary"]

        print(f"\nSequence ID: {result['id']}")
        print(f"Description: {result['description']}")
        print(f"Length: {summary['length']} bp")
        print(f"GC content: {summary['gc_content']}%")
        print(f"QC status: {summary['qc_status']}")

    batch = summarise_batch(results)

    print("\n--- Batch Summary ---")
    print(f"Total records: {batch['total_records']}")
    print(f"Passed: {batch['pass_count']}")
    print(f"Failed: {batch['fail_count']}")
    print(f"Total bases: {batch['total_bases']}")
    print(f"Average length: {batch['average_length']} bp")


if __name__ == "__main__":
    main()
