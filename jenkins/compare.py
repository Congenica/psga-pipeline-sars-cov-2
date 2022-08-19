from pathlib import Path
from typing import Callable, Dict, List, Set
from math import isclose
import pandas as pd

from scripts.util.logging import get_structlog_logger

logger = get_structlog_logger()


class ValidationError(Exception):
    pass


def _check_columns(
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
        # at least one of these two sets will not be empty
        calculated_but_not_expected = sample_names_calc - sample_names_exp
        expected_but_not_calculated = sample_names_exp - sample_names_calc
        logger.error("Found mismatch of samples. In the errors below, focus on these samples only.")
        if calculated_but_not_expected:
            logger.error(f"These unexpected samples were also processed: {calculated_but_not_expected}")
        if expected_but_not_calculated:
            logger.error(f"These expected samples were not calculated: {expected_but_not_calculated}")

    if not columns_to_validate.issubset(df_calc.columns):
        errors = True
        logger.error(f"Columns: {columns_to_validate} are expected in df_calc, but got: {df_calc.columns}")

    if not columns_to_validate.issubset(df_exp.columns):
        errors = True
        logger.error(f"Columns: {columns_to_validate} are expected in df_exp, but got: {df_exp.columns}")

    return errors


def _compare_dataframes(config: Dict, df_calc: pd.DataFrame, df_exp: pd.DataFrame) -> List[str]:
    """
    Compare the tables in the two CSV file paths
    """

    sample_name_column = config["sample_name_column"]
    columns_to_validate = config["columns_to_validate"]
    columns_to_round = config["columns_to_round"]

    errors = _check_columns(df_calc, df_exp, sample_name_column, set(columns_to_validate))

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

    return sample_names


def compare_merged_output_file(
    load_data: Callable, data: Dict, results_csv: Path, expected_results_csv: Path
) -> List[str]:
    """
    Compare the actual vs expected merged output files
    and return the list of samples
    """
    logger.info("Validation of merged output file STARTED")
    df_exp = load_data(data, expected_results_csv)
    df_calc = load_data(data, results_csv)
    sample_names = _compare_dataframes(data, df_calc, df_exp)
    logger.info("Validation PASSED")

    return sample_names


def compare_output_files_set(calc_output_files: Set[Path], exp_output_files: Set[Path]) -> None:
    """
    Compare the expected set of output files against the actual set of output files.
    """
    errors = False

    calc_output_files_only = calc_output_files - exp_output_files
    if calc_output_files_only:
        logger.warning(
            f"These unexpected output files were also generated: {calc_output_files_only}. "
            "IGNORE if the analysis run includes samples which are expected to fail."
        )

    exp_output_files_only = exp_output_files - calc_output_files
    if exp_output_files_only:
        errors = True
        logger.error(f"These expected output files were NOT generated: {exp_output_files_only}")

    if errors:
        raise ValidationError("Validation FAILED. See above for details.")
