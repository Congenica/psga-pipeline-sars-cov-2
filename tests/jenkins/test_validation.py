import pytest
from click.testing import CliRunner

from jenkins.validation import validate_results


@pytest.mark.parametrize(
    "result_file,expected_result_file,tool,exit_code,exception_msg",
    [
        # exact - accept
        (
            "ncov_qc.csv",
            "ncov_qc.csv",
            "ncov2019_artic_nf",
            0,
            None,
        ),
        (
            "all_lineages_report.csv",
            "all_lineages_report.csv",
            "pangolin",
            0,
            None,
        ),
        # approx - accept
        (
            "ncov_qc.csv",
            "ncov_qc_approx.csv",
            "ncov2019_artic_nf",
            0,
            None,
        ),
        (
            "all_lineages_report.csv",
            "all_lineages_report_approx.csv",
            "pangolin",
            0,
            None,
        ),
        # different - reject
        (
            "ncov_qc.csv",
            "ncov_qc_diff.csv",
            "ncov2019_artic_nf",
            1,
            "Validation FAILED. See above for details.",
        ),
        (
            "all_lineages_report.csv",
            "all_lineages_report_diff.csv",
            "pangolin",
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
    tool,
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
            "--tool",
            tool,
        ],
    )

    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert exception_msg in str(rv.exception)
