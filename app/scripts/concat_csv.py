from pathlib import Path
import pandas as pd
import click


def concat(
    input_path: Path,
    output_csv_path: Path,
    sortby_col: str = None,
    ascending: bool = True,
    input_glob_pattern: str = "*.csv",
    input_sep: str = ",",
) -> pd.DataFrame:
    """
    Concatenate the CSV files in input_dir.
    """
    dfs = [pd.read_csv(csv_file, sep=input_sep) for csv_file in input_path.glob(input_glob_pattern)]
    if not dfs:
        raise FileNotFoundError(f"No matching file was found in {input_path}. Cannot concatenate the data frames")
    result = pd.concat(dfs, ignore_index=True)
    if sortby_col:
        if sortby_col not in result.columns:
            raise ValueError(f"Column '{sortby_col}' needed for sorting the data frame not found")
        result.sort_values(by=sortby_col, ascending=ascending, inplace=True)
    result.to_csv(output_csv_path, index=False)
    return result


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
    help="the column to use for sorting the result csv",
)
@click.option(
    "--ascending",
    type=bool,
    default=True,
    help="the direction of the sorting",
)
@click.option(
    "--input-glob-pattern",
    type=str,
    default="*.csv",
    help="The pattern suffix of the input files to concatenate",
)
@click.option(
    "--input-sep",
    type=str,
    default=",",
    help="The separator used in the input files to concatenate",
)
def concat_csv(
    input_path: str,
    output_csv_path: str,
    sortby_col: str,
    ascending: bool,
    input_glob_pattern: str,
    input_sep: str,
) -> None:
    """
    Concatenate the CSV files in the current directory
    """
    concat(Path(input_path), Path(output_csv_path), sortby_col, ascending, input_glob_pattern, input_sep)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    concat_csv()
