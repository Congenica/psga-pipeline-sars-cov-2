import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from scripts.load_pangolin_data_to_db import load_pangolin_data


@pytest.mark.parametrize(
    "sample_name,analysis_run_name,pangolin_sample_columns,pangolin_analysis_run_columns",
    [
        (
            "7174693",
            "just_a_name",
            {
                "pangolin_status": "passed_qc",
                "pangolin_lineage": "AY.98",
                "pangolin_conflict": 0.0,
                "pangolin_ambiguity_score": 1.0,
                "scorpio_call": "Delta (B.1.617.2-like)",
                "scorpio_support": 1.000000,
                "scorpio_conflict": 0.000000,
                "note": "scorpio call: Alt alleles 13; Ref alleles 0; Amb alleles 0; Oth alleles 0",
            },
            {
                "pangolearn_version": "2022-01-20",
                "pangolin_version": "3.1.17",
                "pango_version": "v1.2.123",
            },
        ),
    ],
)
def test_load_pangolin_data_to_db(
    db_session,
    populated_db_session_with_sample,
    test_data_path,
    sample_name,
    analysis_run_name,
    pangolin_sample_columns,
    pangolin_analysis_run_columns,
):
    rv = CliRunner().invoke(
        load_pangolin_data,
        [
            "--pangolin-lineage-report-file",
            test_data_path / f"{sample_name}_lineage_report.csv",
            "--analysis-run-name",
            analysis_run_name,
            "--sample-name",
            sample_name,
        ],
    )

    assert rv.exit_code == 0

    # check pangolin data in sample table
    sample = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            Sample.sample_name == sample_name,
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )
    assert sample is not None
    for col_name, col_val in pangolin_sample_columns.items():
        # get column of sample from string, dynamically
        assert str(getattr(sample, col_name)) == str(col_val)

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
