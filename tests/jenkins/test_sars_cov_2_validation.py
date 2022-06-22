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


def create_output_files(samples: Path, root: Path, ncov_workflow: str, filetype: str):
    df = pd.read_csv(samples)
    sample_names = df["sample_id"].tolist()
    paths = get_expected_output_files(root, sample_names, filetype, ncov_workflow)
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
    "ncov_workflow,filetype",
    [
        (
            "illumina_artic",
            "fastq",
        ),
        (
            "illumina_artic",
            "bam",
        ),
        (
            "medaka_artic",
            "fastq",
        ),
        (
            "no_ncov",
            "fasta",
        ),
    ],
)
@pytest.mark.parametrize(
    "env_var_set,result_file,expected_result_file,exit_code,exception_msg",
    [
        # exact - accept
        (
            True,
            "merged_output.csv",
            "merged_output.csv",
            0,
            None,
        ),
        # approx - accept
        (
            True,
            "merged_output.csv",
            "merged_output_approx.csv",
            0,
            None,
        ),
        # different - reject
        (
            True,
            "merged_output.csv",
            "merged_output_diff.csv",
            1,
            "Validation FAILED. See above for details.",
        ),
        # unset env var
        (
            False,
            "merged_output.csv",
            "merged_output.csv",
            2,
            "Error: Missing option '--psga-output-path'",
        ),
    ],
)
def test_validation(
    monkeypatch,
    tmp_path,
    test_data_path,
    env_var_set,
    result_file,
    expected_result_file,
    ncov_workflow,
    filetype,
    exit_code,
    exception_msg,
):

    if env_var_set:
        monkeypatch.setenv("PSGA_OUTPUT_PATH", str(tmp_path))
    else:
        monkeypatch.delenv("PSGA_OUTPUT_PATH", raising=False)

    actual_path = test_data_path / "jenkins" / "validation" / result_file
    expected_path = test_data_path / "jenkins" / "validation" / expected_result_file

    # here we test the merged output file, but we assume that
    # the output files are as expected
    create_output_files(actual_path, tmp_path, ncov_workflow, filetype)

    rv = CliRunner().invoke(
        validate,
        [
            "--result-path",
            actual_path,
            "--expected-result-path",
            expected_path,
            "sars_cov_2",
            "--ncov-workflow",
            ncov_workflow,
            "--filetype",
            filetype,
        ],
    )
    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert exception_msg in str(rv.exception)

    elif exit_code == 2:
        assert exception_msg in str(rv.output)
