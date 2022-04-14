from pathlib import Path
from typing import Dict, Set, Tuple
import click
import pandas as pd

TOOLS = ["ncov2019_artic_nf", "pangolin"]


class ValidationError(Exception):
    pass


def load_data_frames(
    calculated_csv_path: Path, expected_csv_path: Path, sample_name_column: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load calculated and expected data frames. Perform a basic validation
    """
    df_calc = pd.read_csv(calculated_csv_path)
    df_exp = pd.read_csv(expected_csv_path)

    if sample_name_column not in df_calc.columns or sample_name_column not in df_exp.columns:
        raise ValueError(f"The column {sample_name_column} was not found. Impossible to retrieve samples")

    return df_calc, df_exp


def check_columns(
    df_calc: pd.DataFrame, df_exp: pd.DataFrame, sample_name_column: str, columns_of_interest: Set[str]
) -> Tuple[Set[str], Set[str]]:
    """
    Load sample names and perform basic validations
    """
    sample_names_calc = set(df_calc[sample_name_column])
    sample_names_exp = set(df_exp[sample_name_column])

    if sample_names_calc != sample_names_exp:
        raise ValidationError(f"Expected samples: {sample_names_exp}, but retrieved samples {sample_names_calc}")

    if not columns_of_interest.issubset(df_calc.columns):
        raise ValidationError(f"The columns: {columns_of_interest} are expected in df_calc, but got: {df_calc.columns}")

    if not columns_of_interest.issubset(df_exp.columns):
        raise ValidationError(f"The columns: {columns_of_interest} are expected in df_exp, but got: {df_exp.columns}")

    return sample_names_calc, sample_names_exp


def compare_data_frames(df_calc, df_exp) -> pd.DataFrame:
    """
    Compare two Pandas dataframes and return the diff dataframe
    """
    try:
        return df_calc.compare(df_exp)
    except ValueError as ex:
        raise ValidationError("Cannot compare dataframes due to different labels or shape.") from ex


def validate_csv(config: Dict) -> None:
    """
    Compare the tables in the two CSV file paths
    """
    sample_name_column = config["sample_name_column"]
    coi = config["columns_of_interest"]

    df_calc, df_exp = load_data_frames(config["calculated_output"], config["expected_output"], sample_name_column)
    check_columns(df_calc, df_exp, sample_name_column, set(coi))

    df_calc_sub = df_calc[coi].sort_values(by=[sample_name_column]).round(config["columns_to_round"])
    df_exp_sub = df_exp[coi].sort_values(by=[sample_name_column]).round(config["columns_to_round"])

    df_diff = compare_data_frames(df_calc_sub, df_exp_sub)
    if not df_diff.empty:
        raise ValidationError(f"Calculated results differ from expected results:\n{df_diff}")


validation = {
    # tool which generated the output file to validate
    "ncov2019_artic_nf": {
        # validation function
        "validate": validate_csv,
        # validation function config
        "config": {
            # paths
            "calculated_output": None,
            "expected_output": None,
            # name of the column used for listing the samples
            "sample_name_column": "sample_name",
            # columns to validate
            "columns_of_interest": [
                "sample_name",
                "pct_N_bases",
                "pct_covered_bases",
                "longest_no_N_run",
                "num_aligned_reads",
                "qc_pass",
            ],
            # columns containing floating numbers to round. <colname>: <digits_to_round>
            "columns_to_round": {col: 3 for col in ["pct_N_bases", "pct_covered_bases"]},
        },
    },
    "pangolin": {
        "validate": validate_csv,
        "config": {
            "calculated_output": None,
            "expected_output": None,
            "sample_name_column": "taxon",
            "columns_of_interest": [
                "taxon",
                "lineage",
                "conflict",
                "ambiguity_score",
                "scorpio_call",
                "scorpio_support",
                "scorpio_conflict",
                "scorpio_notes",
                "version",
                "pangolin_version",
                "scorpio_version",
                "constellation_version",
                "is_designated",
                "qc_status",
                "qc_notes",
                "note",
            ],
            "columns_to_round": {
                col: 3 for col in ["conflict", "ambiguity_score", "scorpio_support", "scorpio_conflict"]
            },
        },
    },
}


@click.command()
@click.option(
    "--result-path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The calculated result file",
)
@click.option(
    "--expected-result-path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The expected result file",
)
@click.option(
    "--tool",
    required=True,
    type=click.Choice(TOOLS, case_sensitive=False),
    help="The name of tool generating the data to validate",
)
def validate_results(result_path, expected_result_path, tool):
    """
    Compare the calculated result file against the expected result file.
    """
    validate = validation[tool]["validate"]
    config = validation[tool]["config"]
    config["calculated_output"] = Path(result_path)
    config["expected_output"] = Path(expected_result_path)

    validate(config)

    print("Validation PASSED")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    validate_results()
