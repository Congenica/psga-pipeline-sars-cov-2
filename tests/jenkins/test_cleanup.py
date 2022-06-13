import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from jenkins.cleanup import cleanup, remove_analysis_run


def fetch_sample_names_from_db(db_session, analysis_run_name):
    samples = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .all()
    )
    sample_names = [s.sample_name for s in samples]

    return sample_names


@pytest.mark.parametrize(
    "analysis_run_name,samples",
    [
        (
            "does_not_exist",
            [],
        ),
        (
            "just_a_name",
            ["7284954", "7174693", "8039686", "failed_ncov_qc", "failed_pangolin"],
        ),
    ],
)
def test_remove_analysis_run(populated_db_session_with_sample, db_session, analysis_run_name, samples):
    retrieved_sample_names = fetch_sample_names_from_db(db_session, analysis_run_name)
    assert sorted(samples) == sorted(retrieved_sample_names)

    remove_analysis_run(analysis_run_name)

    retrieved_sample_names = fetch_sample_names_from_db(db_session, analysis_run_name)
    assert [] == sorted(retrieved_sample_names)


@pytest.mark.parametrize(
    "analysis_run_name,samples",
    [
        (
            "just_a_name",
            ["7284954", "7174693", "8039686", "failed_ncov_qc", "failed_pangolin"],
        ),
    ],
)
def test_cleanup(populated_db_session_with_sample, db_session, analysis_run_name, samples):
    retrieved_sample_names = fetch_sample_names_from_db(db_session, analysis_run_name)
    assert sorted(samples) == sorted(retrieved_sample_names)

    rv = CliRunner().invoke(
        cleanup,
        [
            "--analysis-run-name",
            analysis_run_name,
        ],
    )

    assert rv.exit_code == 0

    retrieved_sample_names = fetch_sample_names_from_db(db_session, analysis_run_name)
    assert [] == sorted(retrieved_sample_names)
