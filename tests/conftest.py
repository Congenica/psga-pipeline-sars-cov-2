# pylint: disable=redefined-outer-name
import os
import sys
import tempfile
from pathlib import Path

import pytest
from sqlalchemy import event
from sqlalchemy.orm import sessionmaker

from scripts.db.database import connect

# set this to True to keep changes made to the database while tests are running
# you will be responsible for any cleanup needed before re-running tests in this case
DISABLE_ROLLBACK = False


Session = sessionmaker()


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
    def generator(path, extension, content):
        for sample_id in range(5):
            sample_id = f"NNN{str(sample_id).zfill(5)}"
            fasta_file = path / f"{sample_id}.{extension}"
            fasta_file.write_text(content.format(sample_id=sample_id))

    return generator


@pytest.fixture(autouse=True)
def scripts_test():
    scripts_path = os.path.dirname(__file__).split("/")
    scripts_path[-1] = "scripts"
    sys.path.append("/".join(scripts_path))


@pytest.fixture
def root_genome():
    return Path(__file__).parent.parent / "data" / "FASTA" / "SARS-CoV-2.fasta"


@pytest.fixture
def test_data_path():
    return Path(__file__).parent / "test_data"
