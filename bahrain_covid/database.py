from os import environ

import sqlalchemy
from sqlalchemy.orm import scoped_session, sessionmaker


def connect():
    """
    Connect to the database. Connection details will be retrieved from environment variables.
    """
    host = environ.get("DB_HOST", "localhost")
    port = environ.get("DB_PORT", 5432)

    user = environ.get("DB_USER")
    if not user:
        raise ValueError("DB_USER not set in environment")

    password = environ.get("DB_PASSWORD")
    if not password:
        raise ValueError("DB_PASSWORD not set in environment")

    db_name = environ.get("DB_NAME")
    if not db_name:
        raise ValueError("DB_NAME not set in environment")

    return sqlalchemy.create_engine("postgresql://{user}:{password}@{host}:{port}/{db_name}".format(
        user=user, password=password, host=host, port=port, db_name=db_name
    ))


def session():
    engine = connect()
    session_factory = sessionmaker(expire_on_commit=False)
    Session = scoped_session(session_factory)
    Session.configure(bind=engine)

    return Session
