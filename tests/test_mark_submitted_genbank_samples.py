import pytest
from click.testing import CliRunner

from scripts.db.models import Sample
from scripts.mark_submitted_genbank_samples import mark_submitted_genbank_samples


@pytest.fixture
def populated_db_session_with_samples(db_session):
    sample_names = ["foo", "bar", "buzz"]
    for sample in sample_names:
        db_session.add(Sample(sample_name=sample))

    db_session.commit()
    yield db_session


def test_empty_db(db_session, test_data_path):
    rv = CliRunner().invoke(
        mark_submitted_genbank_samples, ["--sample_names_txt", test_data_path / "samples.txt", "--submit_id", "foo"]
    )

    assert rv.exit_code == 0

    samples = db_session.query(Sample).all()

    assert len(samples) == 0


def test_no_samples(populated_db_session_with_samples, test_data_path):
    db_session = populated_db_session_with_samples
    rv = CliRunner().invoke(
        mark_submitted_genbank_samples, ["--sample_names_txt", test_data_path / "empty.txt", "--submit_id", "foo"]
    )

    assert rv.exit_code == 0

    samples = db_session.query(Sample).filter(Sample.genbank_submit_id.is_(None)).all()
    assert len(samples) == 3


def test_matching_samples(populated_db_session_with_samples, test_data_path):
    db_session = populated_db_session_with_samples
    rv = CliRunner().invoke(
        mark_submitted_genbank_samples, ["--sample_names_txt", test_data_path / "samples.txt", "--submit_id", "foo"]
    )

    assert rv.exit_code == 0

    samples = db_session.query(Sample).filter(Sample.genbank_submit_id == "foo").all()
    assert len(samples) == 3
