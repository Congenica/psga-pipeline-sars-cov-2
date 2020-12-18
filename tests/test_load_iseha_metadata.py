from click.testing import CliRunner

from scripts.db.models import Sample
from scripts.load_iseha_metadata import load_iseha_data


def test_load_good_data(db_session, test_data_path):
    rv = CliRunner().invoke(load_iseha_data, ["--file", test_data_path / "good_iseha_data.tsv"])

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()

    assert len(samples) == 6


def test_load_bad_data(db_session, test_data_path):
    rv = CliRunner().invoke(load_iseha_data, ["--file", test_data_path / "bad_iseha_data.tsv"])

    assert rv.exit_code == 1

    assert (
        rv.output
        == "Invalid row for sample ID HAM44444:\n"
        + "MRN \"HAM111111\" is expected to be an integer, prefixed by 'T' if the person is a foreigner\n"
        + 'AGE "36Y" does not follow any known age formats. It needs to be a '
        "number-only or a number followed by a space then a letter (probably Y)\n"
        + 'GOVERNORATE "MIDDLE" looked right, but not found in database\n'
        + 'AREA "SAKHIR" looked right, but not found in database\n'
        + 'BLOCK "block" should contain digits\n'
        + 'SAMPLE ID "HAM44444" is not a number\n'
        + 'ASSIGN DATE "32/02/2020" is not a valid date: day is out of range for month\n'
        + "Error: Errors encountered: HAM44444\n"
    )

    samples = db_session.query(Sample).all()

    assert len(samples) == 0
