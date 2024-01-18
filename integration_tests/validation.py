from pathlib import Path
import importlib
import click

from integration_tests.loading import load_data_from_csv, get_file_paths
from integration_tests.compare import compare_merged_output_file, compare_output_files_set
from integration_tests.config import data_config
from scripts.check_metadata import SEQUENCING_TECHNOLOGIES
from scripts.util.logger import get_structlog_logger

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)

PATHOGENS = list(data_config)


@click.command()
@click.option(
    "--results-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The calculated result file",
)
@click.option(
    "--expected-results-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The expected result file",
)
@click.option(
    "--output-path",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    required=True,
    help="The PSGA pipeline path where all output files are stored",
)
@click.option(
    "--pathogen",
    required=True,
    type=click.Choice(PATHOGENS, case_sensitive=True),
    help="The pathogen name",
)
@click.option(
    "--sequencing-technology",
    required=False,
    type=click.Choice(SEQUENCING_TECHNOLOGIES, case_sensitive=True),
    default=None,
    help="The name of the sequencing technology",
)
def validate(results_csv: str, expected_results_csv: str, output_path: str, pathogen: str, sequencing_technology: str):
    """
    Compare the calculated result file against the expected result file.
    """

    # Import function based on pathogen module
    # load this lazily as only the module for the invoked pathogen is available in the docker container
    get_expected_output_files = importlib.import_module(f"integration_tests.{pathogen}").get_expected_output_files

    results_csv_path = Path(results_csv)
    expected_results_csv_path = Path(expected_results_csv)
    # do not cast output_path to Path as this can also be "s3://" in get_expected_output_files() below

    if "config" not in data_config[pathogen]:
        raise KeyError(f"key 'config' missing for pathogen '{pathogen}' in data_config")

    validation_config = data_config[pathogen]["config"]

    sample_ids = compare_merged_output_file(
        load_data_from_csv, validation_config, results_csv_path, expected_results_csv_path
    )

    logger.info("Validation of output files set STARTED")
    exp_output_files = get_expected_output_files(output_path, sample_ids, sequencing_technology)
    calc_output_files = get_file_paths(Path(output_path))
    compare_output_files_set(set(calc_output_files), set(exp_output_files))
    logger.info("Validation PASSED")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter,unexpected-keyword-arg
    validate(obj={})
