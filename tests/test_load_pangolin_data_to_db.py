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
            test_data_path / "all_lineages_report.csv",
            "--analysis-run-name",
            analysis_run_name,
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
