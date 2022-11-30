from pathlib import Path
from typing import Dict
import csv
import click


TYPING_SAMPLE_ID_COL = "sample_id"
TYPING_COL = "genotyping"
EXPECTED_TYPING_HEADERS = {
    TYPING_SAMPLE_ID_COL,
    TYPING_COL,
}


def process_typing(
    input_path: Path,
    output_path: Path,
    input_sample_id_col: str,
    input_type_col: str,
) -> None:
    """
    Restructure the information in input_path so that the typing information is stored using
    a dictionary data structure: <sample_id>,{ type1: {...}, type2: {...} }
    """

    with open(input_path, "r", newline="") as input_csv:
        reader = csv.DictReader(input_csv)
        type_keys = [input_sample_id_col, input_type_col]
        output_dict = {}
        typing_dict: Dict[str, Dict] = {}
        for row in reader:
            if not all(key in row for key in type_keys):
                raise KeyError(f"Not all {type_keys} found in {row.keys()}")

            # organise the record values in variables
            sample_id_val = row[input_sample_id_col]
            type_val = row[input_type_col]
            typing_info = {k: v for k, v in row.items() if k not in type_keys}

            if sample_id_val not in output_dict:
                # new sample: reset the typing dict
                typing_dict = {}

            typing_dict[type_val] = typing_info
            output_dict[sample_id_val] = typing_dict

    with open(output_path, "w", newline="") as output_csv:
        writer = csv.DictWriter(output_csv, fieldnames=list(EXPECTED_TYPING_HEADERS))
        writer.writeheader()
        for sample_id, typing_dict in output_dict.items():
            writer.writerow(
                {
                    TYPING_SAMPLE_ID_COL: sample_id,
                    TYPING_COL: typing_dict,
                }
            )


@click.command()
@click.option(
    "--input-path",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="input file containing the sample genotyping",
)
@click.option(
    "--output-csv-path",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the sample id and the typing information as a dictionary",
)
@click.option(
    "--input-sample-id-col",
    type=str,
    required=True,
    help="The sample id",
)
@click.option(
    "--input-type-col",
    type=str,
    required=True,
    help="The column storing the type.",
)
def format_genotyping(
    input_path: str,
    output_csv_path: str,
    input_sample_id_col: str,
    input_type_col: str,
) -> None:
    """
    Structure the typing information associated to the sample using a dictionary
    """
    process_typing(Path(input_path), Path(output_csv_path), input_sample_id_col, input_type_col)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    format_genotyping()
