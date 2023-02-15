from pathlib import Path
import csv
import click


CONTAMINATION_REMOVAL_SAMPLE_ID_COL = "sample_id"
CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL = "contaminated_reads"
EXPECTED_CONTAMINATION_REMOVAL_HEADERS = {
    CONTAMINATION_REMOVAL_SAMPLE_ID_COL,
    CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL,
}


def get_contaminated_reads(input_path: Path) -> int:
    """
    Return the number of reads removed by read-it-and-keep
    """

    def _get_counting(lines: list[str], pattern: str) -> int:
        """e.g. ["Input reads file 1\t75703", "Input reads file 2\t0"] -> 75703"""
        return sum([int(line.split("\t")[1]) for line in lines if line.startswith(pattern)])

    kept_reads = 0
    with open(input_path) as ifd:
        lines = ifd.readlines()
    input_reads = _get_counting(lines, "Input reads")
    kept_reads = _get_counting(lines, "Kept reads")
    contaminated_reads = input_reads - kept_reads
    if contaminated_reads < 0:
        raise ValueError(
            f"Contaminated reads: {contaminated_reads} cannot be negative!"
            f"Input reads: {input_reads}, kept reads: {kept_reads}"
        )

    return contaminated_reads


def write_rik_output_csv(output_csv_path: Path, sample_id: str, contaminated_reads: int) -> None:
    """
    Write a CSV output file for read-it-and-keep
    """
    with open(output_csv_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=list(EXPECTED_CONTAMINATION_REMOVAL_HEADERS))
        writer.writeheader()
        writer.writerow(
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: sample_id,
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: contaminated_reads,
            }
        )


def process_rik(input_path: Path, output_csv_path: Path, sample_id: str) -> None:
    """
    Process the read-it-and-keep output file and store the number of removed reads in a csv file
    """
    contaminated_reads = get_contaminated_reads(input_path)
    write_rik_output_csv(output_csv_path, sample_id, contaminated_reads)


@click.command()
@click.option(
    "--input-path",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="input directory containing the csv files to concatenate",
)
@click.option(
    "--output-csv-path",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file name containing the concatenation of all the csv files in the current directory",
)
@click.option(
    "--sample-id",
    type=str,
    required=True,
    help="The sample id",
)
def contamination_removal(
    input_path: str,
    output_csv_path: str,
    sample_id: str,
) -> None:
    """
    Concatenate the CSV files in the current directory
    """
    process_rik(Path(input_path), Path(output_csv_path), sample_id)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    contamination_removal()
