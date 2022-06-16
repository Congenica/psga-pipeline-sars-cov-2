import pytest
from click.testing import CliRunner

from jenkins.validation import validate_results


@pytest.mark.parametrize(
    "result_file,expected_result_file,pathogen,exit_code,exception_msg",
    [
        # exact - accept
        (
            "merged_output.csv",
            "merged_output.csv",
            "sars_cov_2",
            0,
            None,
        ),
        # approx - accept
        (
            "merged_output.csv",
            "merged_output_approx.csv",
            "sars_cov_2",
            0,
            None,
        ),
        # different - reject
        (
            "merged_output.csv",
            "merged_output_diff.csv",
            "sars_cov_2",
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
    pathogen,
    exit_code,
    exception_msg,
):

    rv = CliRunner().invoke(
        validate_results,
        [
            "--result-path",
            test_data_path / "validation" / result_file,
            "--expected-result-path",
            test_data_path / "validation" / expected_result_file,
            "--pathogen",
            pathogen,
        ],
    )

    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert exception_msg in str(rv.exception)
