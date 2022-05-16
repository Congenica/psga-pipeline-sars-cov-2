import pytest
from click.testing import CliRunner
from pathlib import Path

from scripts.db.models import Sample, SampleQC
from scripts.load_ncov_data_to_db import load_ncov_data
from utils_tests import get_analysis_run_samples


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
                    "qc_pass": True,
                },
                "failed_ncov_qc": {
                    "pct_n_bases": 0.30,
                    "pct_covered_bases": 99.65,
                    "longest_no_n_run": 92725,
                    "num_aligned_reads": 9805873,
                    "qc_pass": False,
                },
                "failed_pangolin": {
                    "pct_n_bases": 0.10,
                    "pct_covered_bases": 95.65,
                    "longest_no_n_run": 21443,
                    "num_aligned_reads": 5851578,
                    "qc_pass": False,
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

    samples_with_no_ncov_qc = Path(tmp_path / "failed_samples.txt")

    rv = CliRunner().invoke(
        load_ncov_data,
        [
            "--ncov-qc-csv-file",
            test_data_path / test_file,
            "--analysis-run-name",
            analysis_run_name,
            "--ncov-qc-depth-directory",
            test_data_path / ncov_qc_depth_dir,
            "--samples-without-qc-file",
            samples_with_no_ncov_qc,
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

    # check failed samples. These are samples where ncov crashed.
    expected_failed_samples = [k for k, v in expected_values.items() if not v["qc_pass"]]
    with open(samples_with_no_ncov_qc, "r") as ifr:
        failed_samples = ifr.read().splitlines()
    assert sorted(expected_failed_samples) == sorted(failed_samples)
