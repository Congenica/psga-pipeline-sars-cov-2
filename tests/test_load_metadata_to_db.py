from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from scripts.load_metadata_to_db import load_metadata


def test_load_good_data(db_session, test_data_path):
    analysis_name = "success"
    rv = CliRunner().invoke(
        load_metadata,
        [
            "--file",
            test_data_path / "good_metadata.tsv",
            "--analysis-run-name",
            analysis_name,
            "--pipeline-version",
            "1.0.0",
        ],
    )

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()
    assert len(samples) == 6

    analysis_run = db_session.query(AnalysisRun).one_or_none()
    assert analysis_run.analysis_run_name == analysis_name


def test_load_bad_data(db_session, test_data_path):
    analysis_name = "unsuccess"
    rv = CliRunner().invoke(
        load_metadata,
        [
            "--file",
            test_data_path / "bad_metadata.tsv",
            "--analysis-run-name",
            analysis_name,
            "--pipeline-version",
            "1.0.0",
        ],
    )

    assert rv.exit_code == 1

    assert (
        rv.output
        == "Invalid row for sample ID HAM44444:\n"
        + 'ASSIGN DATE "32/02/2020" is not a valid date: day is out of range for month\n'
        + "Error: Errors encountered: HAM44444\n"
    )

    samples = db_session.query(Sample).all()
    assert len(samples) == 0

    analysis_run = db_session.query(AnalysisRun).one_or_none()
    assert analysis_run.analysis_run_name == analysis_name
