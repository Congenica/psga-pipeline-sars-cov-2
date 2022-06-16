from pathlib import Path
from typing import Dict, Set
from math import isclose
import click
import pandas as pd

from scripts.util.logging import get_structlog_logger
from jenkins.config import data_config

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)


class ValidationError(Exception):
    pass


def check_columns(
    df_calc: pd.DataFrame, df_exp: pd.DataFrame, sample_name_column: str, columns_to_validate: Set[str]
) -> bool:
    """
    Load sample names and perform basic validations
    """
    sample_names_calc = set(df_calc[sample_name_column])
    sample_names_exp = set(df_exp[sample_name_column])
    errors = False

    if sample_names_calc != sample_names_exp:
        errors = True
        logger.error(f"Expected samples: {sample_names_exp}, but retrieved samples {sample_names_calc}")

    if not columns_to_validate.issubset(df_calc.columns):
        errors = True
        logger.error(f"Columns: {columns_to_validate} are expected in df_calc, but got: {df_calc.columns}")

    if not columns_to_validate.issubset(df_exp.columns):
        errors = True
        logger.error(f"Columns: {columns_to_validate} are expected in df_exp, but got: {df_exp.columns}")

    return errors


def validate(config: Dict, df_calc: pd.DataFrame, df_exp: pd.DataFrame) -> None:
    """
    Compare the tables in the two CSV file paths
    """

    sample_name_column = config["sample_name_column"]
    columns_to_validate = config["columns_to_validate"]
    columns_to_round = config["columns_to_round"]

    errors = check_columns(df_calc, df_exp, sample_name_column, set(columns_to_validate))

    sample_names = df_calc[sample_name_column].tolist()

    for col_to_compare in columns_to_validate:
        col_calc = df_calc[col_to_compare].fillna(0).tolist()
        col_exp = df_exp[col_to_compare].fillna(0).tolist()

        col_to_round = col_to_compare in columns_to_round
        diff_data = {}

        for sn, calc, exp in zip(sample_names, col_calc, col_exp):
            # use math.isclose instead of numpy.isclose because the former is symmetric in a and b ; e.g.
            # abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol) , instead of
            # abs(a-b) <= (atol + rtol * abs(b))
            # Our expected dataframe is from an execution. The values are not averages
            match = (
                isclose(calc, exp, abs_tol=columns_to_round[col_to_compare])
                if col_to_compare in columns_to_round
                else calc == exp
            )
            if not match:
                diff_data[sn] = {"calculated": calc, "expected": exp}

        if diff_data:
            errors = True
            abstol_info = f" using abs_tol {columns_to_round[col_to_compare]}" if col_to_round else ""
            logger.error(f"Found mismatches for column: {col_to_compare}. Calculated vs Expected values{abstol_info}:")
            for sample, comp in diff_data.items():
                logger.error(sample=sample, calculated=comp["calculated"], expected=comp["expected"])

    if errors:
        raise ValidationError("Validation FAILED. See above for details.")


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
    "--pathogen",
    required=True,
    type=click.Choice(data_config, case_sensitive=False),
    help="The name of the pathogen",
)
def validate_results(result_path, expected_result_path, pathogen):
    """
    Compare the calculated result file against the expected result file.
    """
    data = data_config[pathogen]["config"]

    logger.info("Validation of output file STARTED")
    load_data = data_config[pathogen]["load_data"]
    df_exp = load_data(data, Path(expected_result_path))
    df_calc = load_data(data, Path(result_path))
    validate(data, df_calc, df_exp)
    logger.info("Validation PASSED")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    validate_results()
