from pathlib import Path
import pytest
from click.testing import CliRunner
import pandas as pd

from jenkins.compare import compare_merged_output_file, ValidationError
from jenkins.config import data_config
from jenkins.loading import load_data_from_csv
from jenkins.sars_cov_2 import get_expected_output_files
from jenkins.validation import validate

from tests.jenkins.util import create_paths


def create_output_files(samples: Path, root: Path, sequencing_technology: str):
    df = pd.read_csv(samples)
    sample_names = df["sample_id"].tolist()
    paths = get_expected_output_files(root, sample_names, sequencing_technology)
    create_paths(paths)


@pytest.mark.parametrize(
    "result_file,expected_result_file,exc",
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
    test_data_path,
    result_file,
    expected_result_file,
    exc,
):
    actual_path = test_data_path / "jenkins" / "validation" / result_file
    expected_path = test_data_path / "jenkins" / "validation" / expected_result_file
    data = data_config["sars_cov_2"]["config"]

    if exc:
        with pytest.raises(ValidationError, match=exc):
            compare_merged_output_file(load_data_from_csv, data, actual_path, expected_path)
    else:
        sample_names = compare_merged_output_file(load_data_from_csv, data, actual_path, expected_path)
        df = pd.read_csv(expected_path)
        expected_sample_names = df["sample_id"].tolist()
        assert set(sample_names) == set(expected_sample_names)


@pytest.mark.parametrize(
    "sequencing_technology",
    ["illumina", "ont", "unknown"],
)
@pytest.mark.parametrize(
    "result_file,expected_result_file,exit_code,exception_msg",
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
    test_data_path,
    result_file,
    expected_result_file,
    sequencing_technology,
    exit_code,
    exception_msg,
):

    actual_path = test_data_path / "jenkins" / "validation" / result_file
    expected_path = test_data_path / "jenkins" / "validation" / expected_result_file

    # here we test the merged output file, but we assume that
    # the output files are as expected
    create_output_files(actual_path, tmp_path, sequencing_technology)

    rv = CliRunner().invoke(
        validate,
        [
            "--result-path",
            actual_path,
            "--expected-result-path",
            expected_path,
            "--output-path",
            tmp_path,
            "sars_cov_2",
            "--sequencing-technology",
            sequencing_technology,
        ],
    )
    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert exception_msg in str(rv.exception)

    elif exit_code == 2:
        assert exception_msg in str(rv.output)
