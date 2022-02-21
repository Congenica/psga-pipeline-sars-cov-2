import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample, SampleQC
from scripts.load_ncov_data_to_db import load_ncov_data


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
            },
        ),
    ],
)
def test_load_pangolin_data_to_db(
    db_session,
    populated_db_session_with_sample,
    test_data_path,
    test_file,
    analysis_run_name,
    ncov_qc_depth_dir,
    expected_values,
):
    rv = CliRunner().invoke(
        load_ncov_data,
        [
            "--ncov-qc-csv-file",
            test_data_path / test_file,
            "--analysis-run-name",
            analysis_run_name,
            "--ncov-qc-depth-directory",
            test_data_path / ncov_qc_depth_dir,
        ],
    )

    assert rv.exit_code == 0

    # check pangolin data in sample table
    samples = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .all()
    )
    assert len(samples) == 3
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
        assert sample_qc is not None
        for col_name, col_val in expected_values[sample_name].items():
            # get column of sample from string, dynamically
            assert getattr(sample_qc, col_name) == col_val
