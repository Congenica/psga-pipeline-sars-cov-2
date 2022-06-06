from pathlib import Path
import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from scripts.sars_cov_2.load_ncov_results_to_db import load_ncov_results_to_db
from scripts.sars_cov_2.load_pangolin_results_to_db import load_pangolin_results_to_db
from jenkins.validation import validate_results


@pytest.fixture
def validation_analysis_run(db_session):
    analysis_run_name = "validation_analysis_run"
    db_session.add(AnalysisRun(analysis_run_name=analysis_run_name))
    analysis_run = (
        db_session.query(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )

    for sn in [
        "ERR5069314",
        "ERR5069318",
        "ERR5069325",
        "ERR5440057",
        "ERR5440118",
        "ERR5440120",
        "ERR5440121",
        "ERR5440122",
        "ERR5440123",
        "ERR5465776",
        "ERR5465779",
        "ERR5465780",
        "ERR5465785",
        "ERR5465787",
        "ERR5465789",
        "ERR5466036",
        "ERR5466039",
        "ERR5466041",
        "ERR5466043",
        "ERR5853949",
        "ERR8269258",
        "ERR8269259",
        "ERR8269260",
        "ERR8269268",
        "ERR8269269",
        "ERR8269270",
        "ERR8269276",
        "ERR8269277",
        "ERR8269541",
        "ERR8269542",
    ]:
        db_session.add(
            Sample(
                sample_name=sn,
                analysis_run_id=analysis_run.analysis_run_id,
                metadata_loaded=True,
            )
        )

    db_session.commit()
    yield db_session


@pytest.mark.parametrize(
    "result_file,expected_result_file,tool,analysis_run_name,exit_code,exception_msg",
    [
        # exact - accept
        (
            "ncov_qc.csv",
            "ncov_qc.csv",
            "ncov2019_artic_nf",
            "validation_analysis_run",
            0,
            None,
        ),
        (
            "all_lineages_report.csv",
            "all_lineages_report.csv",
            "pangolin",
            "validation_analysis_run",
            0,
            None,
        ),
        # approx - accept
        (
            "ncov_qc.csv",
            "ncov_qc_approx.csv",
            "ncov2019_artic_nf",
            "validation_analysis_run",
            0,
            None,
        ),
        (
            "all_lineages_report.csv",
            "all_lineages_report_approx.csv",
            "pangolin",
            "validation_analysis_run",
            0,
            None,
        ),
        # different - reject
        (
            "ncov_qc.csv",
            "ncov_qc_diff.csv",
            "ncov2019_artic_nf",
            "validation_analysis_run",
            1,
            "Validation FAILED. See above for details.",
        ),
        (
            "all_lineages_report.csv",
            "all_lineages_report_diff.csv",
            "pangolin",
            "validation_analysis_run",
            1,
            "Validation FAILED. See above for details.",
        ),
    ],
)
def test_validation(
    tmp_path,
    test_data_path,
    validation_analysis_run,
    result_file,
    expected_result_file,
    tool,
    analysis_run_name,
    exit_code,
    exception_msg,
):

    # load the result_file to the DB so that we can test this too
    if tool == "ncov2019_artic_nf":
        rv = CliRunner().invoke(
            load_ncov_results_to_db,
            [
                "--ncov-qc-csv-file",
                test_data_path / "validation" / result_file,
                "--analysis-run-name",
                analysis_run_name,
                "--ncov-qc-depth-directory",
                test_data_path / "validation" / "qc_plots",
                "--samples-without-ncov-qc-file",
                Path(tmp_path / "no_ncov_qc_samples.txt"),
                "--samples-with-failed-ncov-qc-file",
                Path(tmp_path / "failed_ncov_qc_samples.txt"),
                "--samples-with-passed-ncov-qc-file",
                Path(tmp_path / "passed_ncov_qc_samples.txt"),
            ],
        )
    elif tool == "pangolin":
        rv = CliRunner().invoke(
            load_pangolin_results_to_db,
            [
                "--pangolin-lineage-report-file",
                test_data_path / "validation" / result_file,
                "--analysis-run-name",
                analysis_run_name,
                "--samples-with-unknown-pangolin-status-file",
                Path(tmp_path / "unknown_pangolin.txt"),
                "--samples-with-failed-pangolin-status-file",
                Path(tmp_path / "failed_pangolin.txt"),
                "--samples-with-passed-pangolin-status-file",
                Path(tmp_path / "passed_pangolin.txt"),
            ],
        )
    else:
        raise ValueError("<tool> was not recognised")

    assert rv.exit_code == 0

    rv = CliRunner().invoke(
        validate_results,
        [
            "--result-path",
            test_data_path / "validation" / result_file,
            "--expected-result-path",
            test_data_path / "validation" / expected_result_file,
            "--tool",
            tool,
            "--analysis-run-name",
            analysis_run_name,
        ],
    )

    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert exception_msg in str(rv.exception)
