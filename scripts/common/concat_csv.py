from pathlib import Path
import pandas as pd
import click


def concat(input_path: Path, output_csv_path: Path, sortby_col: str) -> None:
    """
    Concatenate the CSV files in input_dir.
    """
    dfs = [pd.read_csv(csv_file) for csv_file in input_path.glob("*.csv")]
    result = pd.concat(dfs, ignore_index=True)
    if sortby_col not in result.columns:
        raise ValueError(f"Column '{sortby_col}' needed for sorting the data frame was not found")
    result.sort_values(by=[sortby_col], inplace=True)
    result.to_csv(output_csv_path, index=False)


@click.command()
@click.option(
    "--input-path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
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
    "--sortby-col",
    type=str,
    required=True,
    help="the name of the column to use for sorting the result csv",
)
def concat_csv(
    input_path: str,
    output_csv_path: str,
    sortby_col: str,
) -> None:
    """
    Concatenate the CSV files in the current directory
    """
    concat(Path(input_path), Path(output_csv_path), sortby_col)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    concat_csv()
