import os
from pathlib import Path

from click.testing import CliRunner

from scripts.db.models import Sample
from scripts.load_iseha_metadata import load_iseha_data

TEST_DATA_PATH = Path(os.path.join(os.path.dirname(__file__), "test_data"))


def test_load_good_data(db_session):
    rv = CliRunner().invoke(load_iseha_data, ["--file", open(TEST_DATA_PATH / "good_iseha_data.tsv", "r")])

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()

    assert len(samples) == 4


def test_load_bad_data(db_session):
    rv = CliRunner().invoke(load_iseha_data, ["--file", open(TEST_DATA_PATH / "bad_iseha_data.tsv", "r")])

    assert rv.exit_code == 0

    assert (
        rv.output
        == "Invalid row for sample ID HAM44444:\n"
        + 'MRN "HAM111111" is not an integer\n'
        + 'AGE "36" should be a number followed by a space then a letter (probably Y)\n'
        + 'GOVERNORATE "MIDDLE" looked right, but not found in database\n'
        + 'AREA "SAKHIR" looked right, but not found in database\n'
        + 'BLOCK "block" should contain digits\n'
        + 'SAMPLE ID "HAM44444" is not a number\n'
        + 'ASSIGN DATE "32/02/2020" is not a valid date: day is out of range for month\n'
        + "Inserted samples: \n"
        + "Updated samples: \n"
        + "Errors encountered, these samples were not inserted or updated: HAM44444\n"
    )

    samples = db_session.query(Sample).all()

    assert len(samples) == 0
