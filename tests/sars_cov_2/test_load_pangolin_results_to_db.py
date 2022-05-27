from pathlib import Path
import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun
from scripts.sars_cov_2.load_pangolin_results_to_db import load_pangolin_results_to_db
from utils_tests import get_analysis_run_samples, read_samples_from_file


def check_pangolin_status(input_path, samples, status):
    processed_samples = read_samples_from_file(input_path)
    expected_samples = [k for k, v in samples.items() if v["pangolin_status"] == status]
    assert sorted(expected_samples) == sorted(processed_samples)


@pytest.mark.parametrize(
    "test_file,analysis_run_name,pangolin_sample_columns,pangolin_analysis_run_columns",
    [
        (
            "all_lineages_report.csv",
            "just_a_name",
            {
                "7174693": {
                    "sample_name": "7174693",
                    "pangolin_status": "PASS",
                    "pangolin_lineage": "AY.98",
                    "pangolin_conflict": 0.0,
                    "pangolin_ambiguity_score": None,
                    "scorpio_call": "Delta (B.1.617.2-like)",
                    "scorpio_support": 1.000000,
                    "scorpio_conflict": 0.000000,
                    "scorpio_notes": "scorpio call: Alt alleles 13; Ref alleles 0; Amb alleles 0; Oth alleles 0",
                    "is_designated": False,
                    "qc_notes": "Ambiguous_content:0.08",
                    "note": "Usher placements: AY.98(1/1)",
                },
                "8039686": {
                    "sample_name": "8039686",
                    "pangolin_status": "FAIL",
                },
                "7284954": {
                    "sample_name": "7284954",
                    "pangolin_status": "FAIL",
                },
                "failed_ncov_qc": {
                    "sample_name": "failed_ncov_qc",
                    "pangolin_status": "UNKNOWN",
                },
                "failed_pangolin": {
                    "sample_name": "failed_pangolin",
                    "pangolin_status": "UNKNOWN",
                },
            },
            {
                "pangolin_version": "4.0.2",
                "pangolin_data_version": "PUSHER-v1.2.133",
                "constellation_version": "v0.1.4",
                "scorpio_version": "0.3.16",
            },
        ),
    ],
)
def test_load_pangolin_data_to_db(
    db_session,
    populated_db_session_with_sample,
    tmp_path,
    test_data_path,
    test_file,
    analysis_run_name,
    pangolin_sample_columns,
    pangolin_analysis_run_columns,
):

    samples_with_unknown_pangolin_status = Path(tmp_path / "unknown_pangolin.txt")
    samples_with_failed_pangolin_status = Path(tmp_path / "failed_pangolin.txt")
    samples_with_passed_pangolin_status = Path(tmp_path / "passed_pangolin.txt")

    rv = CliRunner().invoke(
        load_pangolin_results_to_db,
        [
            "--pangolin-lineage-report-file",
            test_data_path / test_file,
            "--analysis-run-name",
            analysis_run_name,
            "--samples-with-unknown-pangolin-status-file",
            samples_with_unknown_pangolin_status,
            "--samples-with-failed-pangolin-status-file",
            samples_with_failed_pangolin_status,
            "--samples-with-passed-pangolin-status-file",
            samples_with_passed_pangolin_status,
        ],
    )

    assert rv.exit_code == 0

    # check samples
    samples = get_analysis_run_samples(db_session, analysis_run_name)
    for sample in samples:
        sample_name = sample.sample_name
        if sample_name in pangolin_sample_columns:
            # samples with pangolin status: fail or pass
            for col_name, col_val in pangolin_sample_columns[sample_name].items():
                # get column of sample from string, dynamically
                assert str(getattr(sample, col_name)) == str(col_val)

    # check notification files
    check_pangolin_status(samples_with_unknown_pangolin_status, pangolin_sample_columns, "UNKNOWN")
    check_pangolin_status(samples_with_failed_pangolin_status, pangolin_sample_columns, "FAIL")
    check_pangolin_status(samples_with_passed_pangolin_status, pangolin_sample_columns, "PASS")

    # check pangolin data in analysis_run table
    analysis_run = (
        db_session.query(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )
    assert analysis_run is not None
    for col_name, col_val in pangolin_analysis_run_columns.items():
        # get column of sample from string, dynamically
        assert getattr(analysis_run, col_name) == col_val
