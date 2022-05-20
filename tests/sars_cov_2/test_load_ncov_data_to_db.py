import pytest
from click.testing import CliRunner
from pathlib import Path

from scripts.db.models import Sample, SampleQC
from scripts.sars_cov_2.load_ncov_data_to_db import load_ncov_data
from utils_tests import get_analysis_run_samples, read_samples_from_file


def check_ncov_qc(input_path, samples, status):
    processed_samples = read_samples_from_file(input_path)
    expected_samples = [k for k, v in samples.items() if v["qc_pass"] == status]
    assert sorted(expected_samples) == sorted(processed_samples)


@pytest.mark.parametrize(
    "test_file,analysis_run_name,ncov_qc_depth_dir,expected_values",
    [
        (
            "ncov_test.qc.csv",
            "just_a_name",
            "qc_plots",
            {
                "7284954": {
                    "pct_n_bases": 0.30,
                    "pct_covered_bases": 99.67,
                    "longest_no_n_run": 29798,
                    "num_aligned_reads": 2406868,
                    "qc_pass": True,
                },
                "7174693": {
                    "pct_n_bases": 0.12,
                    "pct_covered_bases": 99.86,
                    "longest_no_n_run": 29822,
                    "num_aligned_reads": 2989521,
                    "qc_pass": True,
                },
                "8039686": {
                    "pct_n_bases": 0.30,
                    "pct_covered_bases": 99.65,
                    "longest_no_n_run": 29783,
                    "num_aligned_reads": 3805878,
                    "qc_pass": False,
                },
                "failed_ncov_qc": {
                    "qc_pass": None,
                },
                "failed_pangolin": {
                    "qc_pass": None,
                },
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
    ncov_qc_depth_dir,
    expected_values,
):

    samples_with_no_ncov_qc = Path(tmp_path / "no_ncov_qc_samples.txt")
    samples_with_failed_ncov_qc = Path(tmp_path / "failed_ncov_qc_samples.txt")
    samples_with_passed_ncov_qc = Path(tmp_path / "passed_ncov_qc_samples.txt")

    rv = CliRunner().invoke(
        load_ncov_data,
        [
            "--ncov-qc-csv-file",
            test_data_path / test_file,
            "--analysis-run-name",
            analysis_run_name,
            "--ncov-qc-depth-directory",
            test_data_path / ncov_qc_depth_dir,
            "--samples-without-ncov-qc-file",
            samples_with_no_ncov_qc,
            "--samples-with-failed-ncov-qc-file",
            samples_with_failed_ncov_qc,
            "--samples-with-passed-ncov-qc-file",
            samples_with_passed_ncov_qc,
        ],
    )

    assert rv.exit_code == 0

    samples = get_analysis_run_samples(db_session, analysis_run_name)
    for sample in samples:
        sample_name = sample.sample_name
        sample_qc = (
            db_session.query(SampleQC)
            .join(Sample)
            .filter(
                Sample.sample_name == sample.sample_name,
            )
            .one_or_none()
        )
        # check successfull samples
        if sample_qc is not None:
            for col_name, col_val in expected_values[sample_name].items():
                # get column of sample from string, dynamically
                assert getattr(sample_qc, col_name) == col_val

    check_ncov_qc(samples_with_no_ncov_qc, expected_values, None)
    check_ncov_qc(samples_with_failed_ncov_qc, expected_values, False)
    check_ncov_qc(samples_with_passed_ncov_qc, expected_values, True)
