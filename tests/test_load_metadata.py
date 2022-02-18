from click.testing import CliRunner

from scripts.db.models import Sample
from scripts.load_metadata import load_metadata


def test_load_good_data(db_session, test_data_path):
    rv = CliRunner().invoke(load_metadata, ["--file", test_data_path / "good_metadata.tsv"])

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()

    assert len(samples) == 6


def test_load_bad_data(db_session, test_data_path):
    rv = CliRunner().invoke(load_metadata, ["--file", test_data_path / "bad_metadata.tsv"])

    assert rv.exit_code == 1

    assert (
        rv.output
        == "Invalid row for sample ID HAM44444:\n"
        + 'ASSIGN DATE "32/02/2020" is not a valid date: day is out of range for month\n'
        + "Error: Errors encountered: HAM44444\n"
    )

    samples = db_session.query(Sample).all()

    assert len(samples) == 0
