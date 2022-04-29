# pylint: disable=redefined-outer-name
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import pytest
from pytest_socket import disable_socket
from scripts.db.database import connect
from scripts.db.models import AnalysisRun, Sample
from sqlalchemy import event
from sqlalchemy.orm import sessionmaker

# set this to True to keep changes made to the database while tests are running
# you will be responsible for any cleanup needed before re-running tests in this case
DISABLE_ROLLBACK = False


Session = sessionmaker()


def pytest_runtest_setup():
    disable_socket()


@pytest.fixture
def db_engine():
    return connect()


@pytest.fixture
def db_session(monkeypatch, db_engine):
    connection = db_engine.connect()
    tx = connection.begin()

    session = Session(bind=connection)

    if not DISABLE_ROLLBACK:
        session.begin_nested()

    @event.listens_for(session, "after_transaction_end")
    def restart_savepoint(sess, transaction):
        if transaction.nested and not transaction._parent.nested and not DISABLE_ROLLBACK:
            sess.expire_all()
            sess.begin_nested()

    monkeypatch.setattr("scripts.db.database.create_session", lambda: session)

    yield session

    if not DISABLE_ROLLBACK:
        session.close()
        tx.rollback()
    else:
        tx.commit()
        session.commit()

    connection.close()


@pytest.fixture
def tmp_path():
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def fasta_file_generator():
    # generate consensus fasta files
    def generator(path, extension, content):
        for sample_id in ["1", "ERR12313", "ERR12313_barcode001", "529d82ab-51b7-433c-8be6-f1d59759b943"]:
            sample_id = f"NNN{str(sample_id).zfill(5)}"
            fasta_file = path / f"{sample_id}.consensus.{extension}"
            fasta_file.write_text(content.format(sample_id=sample_id))

    return generator


@pytest.fixture(autouse=True)
def scripts_test():
    scripts_path = os.path.dirname(__file__).split("/")
    scripts_path[-1] = "scripts"
    sys.path.append("/".join(scripts_path))


@pytest.fixture
def test_data_path():
    return Path(__file__).parent / "test_data"


@pytest.fixture
def test_data_path_genbank_input(test_data_path):
    return test_data_path / "genbank_input"


@pytest.fixture
def test_data_path_genbank_reference(test_data_path):
    return test_data_path / "genbank_reference"


@pytest.fixture
def db_fetcher_by_name(db_session):
    def fetcher(db_object, name):
        return db_session.query(db_object).filter(db_object.name == name).one_or_none()

    return fetcher


@pytest.fixture
def sample_generator(db_session, db_fetcher_by_name):
    def generate_sample(lineage="B.1", date_collected=datetime.now()):
        return db_session.add(
            Sample(
                date_collected=date_collected,
                pangolin_lineage=lineage,
            )
        )

    return generate_sample


@pytest.fixture
def populated_db_session_with_sample(db_session):
    analysis_run_name = "just_a_name"
    db_session.add(AnalysisRun(analysis_run_name=analysis_run_name))
    analysis_run = (
        db_session.query(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )

    for sn in ["7284954", "7174693", "8039686"]:
        db_session.add(
            Sample(
                sample_name=sn,
                analysis_run_id=analysis_run.analysis_run_id,
                metadata_loaded=True,
            )
        )

    db_session.commit()
    yield db_session
