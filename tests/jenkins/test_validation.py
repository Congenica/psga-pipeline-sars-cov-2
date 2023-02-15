from pathlib import Path
import importlib
import pytest
from click.testing import CliRunner
import pandas as pd

from scripts.util.metadata import SAMPLE_ID
from jenkins.compare import compare_merged_output_file, ValidationError
from jenkins.config import data_config
from jenkins.loading import load_data_from_csv
from jenkins.validation import validate

from tests.jenkins.util import create_paths


def create_output_files(pathogen: str, samples: Path, root: Path, sequencing_technology: str):
    df = pd.read_csv(samples)
    sample_names = df[SAMPLE_ID].tolist()
    # Import function based on pathogen module
    # load this lazily as only the module for the invoked pathogen is available in the docker container
    get_expected_output_files = importlib.import_module(f"jenkins.{pathogen}").get_expected_output_files
    paths = get_expected_output_files(str(root), sample_names, sequencing_technology)
    create_paths(paths)


@pytest.mark.parametrize("pathogen", ["sars_cov_2", "synthetic"])
@pytest.mark.parametrize(
    "results_csv,expected_results_csv,exc",
    [
        # exact - accept
        (
            "merged_output.csv",
            "merged_output.csv",
            None,
        ),
        # approx - accept
        (
            "merged_output.csv",
            "merged_output_approx.csv",
            None,
        ),
        # different - reject
        (
            "merged_output.csv",
            "merged_output_diff.csv",
            "Validation FAILED. See above for details.",
        ),
    ],
)
def test_compare_merged_output_file(
    integration_test_validation_data_path,
    pathogen,
    results_csv,
    expected_results_csv,
    exc,
):
    actual_path = integration_test_validation_data_path / pathogen / results_csv
    expected_path = integration_test_validation_data_path / pathogen / expected_results_csv
    data = data_config[pathogen]["config"]

    if exc:
        with pytest.raises(ValidationError, match=exc):
            compare_merged_output_file(load_data_from_csv, data, actual_path, expected_path)
    else:
        sample_names = compare_merged_output_file(load_data_from_csv, data, actual_path, expected_path)
        df = pd.read_csv(expected_path)
        expected_sample_names = df[SAMPLE_ID].tolist()
        assert set(sample_names) == set(expected_sample_names)


@pytest.mark.parametrize("pathogen", ["sars_cov_2", "synthetic"])
@pytest.mark.parametrize(
    "sequencing_technology",
    ["illumina", "ont", "unknown"],
)
@pytest.mark.parametrize(
    "results_csv,expected_results_csv,exit_code,exception_msg",
    [
        # exact - accept
        (
            "merged_output.csv",
            "merged_output.csv",
            0,
            None,
        ),
        # approx - accept
        (
            "merged_output.csv",
            "merged_output_approx.csv",
            0,
            None,
        ),
        # different - reject
        (
            "merged_output.csv",
            "merged_output_diff.csv",
            1,
            "Validation FAILED. See above for details.",
        ),
    ],
)
def test_validation(
    tmp_path,
    integration_test_validation_data_path,
    pathogen,
    results_csv,
    expected_results_csv,
    sequencing_technology,
    exit_code,
    exception_msg,
):

    actual_path = integration_test_validation_data_path / pathogen / results_csv
    expected_path = integration_test_validation_data_path / pathogen / expected_results_csv

    # here we test the merged output file, but we assume that
    # the output files are as expected
    create_output_files(pathogen, actual_path, tmp_path, sequencing_technology)

    rv = CliRunner().invoke(
        validate,
        [
            "--results-csv",
            actual_path,
            "--expected-results-csv",
            expected_path,
            "--output-path",
            tmp_path,
            "--pathogen",
            pathogen,
            "--sequencing-technology",
            sequencing_technology,
        ],
    )
    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert exception_msg in str(rv.exception)

    elif exit_code == 2:
        assert exception_msg in str(rv.output)
