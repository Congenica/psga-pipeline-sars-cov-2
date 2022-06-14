import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from scripts.genbank.mark_submitted_genbank_samples import mark_submitted_genbank_samples

ANALYSIS_RUN_NAME = "just_a_name"


@pytest.fixture
def populated_db_session_with_samples(db_session):
    db_session.add(AnalysisRun(analysis_run_name=ANALYSIS_RUN_NAME))
    analysis_run = (
        db_session.query(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == ANALYSIS_RUN_NAME,
        )
        .one_or_none()
    )

    sample_names = ["foo", "bar", "buzz"]
    for sample in sample_names:
        db_session.add(Sample(analysis_run_id=analysis_run.analysis_run_id, sample_name=sample))

    db_session.commit()
    yield db_session


def test_empty_db(db_session, test_data_path):
    rv = CliRunner().invoke(
        mark_submitted_genbank_samples,
        [
            "--analysis-run-name",
            ANALYSIS_RUN_NAME,
            "--sample-names-txt",
            test_data_path / "samples.txt",
            "--submit-id",
            "foo",
        ],
    )

    assert rv.exit_code == 0

    samples = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == ANALYSIS_RUN_NAME,
        )
        .all()
    )

    assert len(samples) == 0


def test_no_samples(populated_db_session_with_samples, test_data_path):
    db_session = populated_db_session_with_samples
    rv = CliRunner().invoke(
        mark_submitted_genbank_samples,
        [
            "--analysis-run-name",
            ANALYSIS_RUN_NAME,
            "--sample-names-txt",
            test_data_path / "empty.txt",
            "--submit-id",
            "foo",
        ],
    )

    assert rv.exit_code == 0

    samples = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == ANALYSIS_RUN_NAME,
            Sample.genbank_submit_id.is_(None),
        )
        .all()
    )
    assert len(samples) == 3


def test_matching_samples(populated_db_session_with_samples, test_data_path):
    db_session = populated_db_session_with_samples
    rv = CliRunner().invoke(
        mark_submitted_genbank_samples,
        [
            "--analysis-run-name",
            ANALYSIS_RUN_NAME,
            "--sample-names-txt",
            test_data_path / "samples.txt",
            "--submit-id",
            "foo",
        ],
    )

    assert rv.exit_code == 0

    samples = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == ANALYSIS_RUN_NAME,
            Sample.genbank_submit_id == "foo",
        )
        .all()
    )
    assert len(samples) == 3
