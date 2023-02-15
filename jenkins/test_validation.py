from click.testing import CliRunner

from jenkins.validation import validate

# call:
# pytest test_validation.py --results-csv results.csv --expected-results-csv results2.csv \
#        --output-path /tmp --pathogen sars_cov_2 --sequencing-technology ont


def test_validation(
    results_csv: str, expected_results_csv: str, pathogen: str, output_path: str, sequencing_technology: str
):
    """
    Run the validation script via pytest
    """

    cmd_config = [
        "--results-csv",
        results_csv,
        "--expected-results-csv",
        expected_results_csv,
        "--output-path",
        output_path,
        "--pathogen",
        pathogen,
        "--sequencing-technology",
        sequencing_technology,
    ]

    rv = CliRunner().invoke(
        validate,
        cmd_config,
    )

    # leave this so that we can get the detailed output of the test
    print(rv.output)

    assert rv.exit_code == 0
