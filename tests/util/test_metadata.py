import uuid

from scripts.util.metadata import (
    is_valid_uuid,
)


def test_is_valid_uuid():
    uuid_str = str(uuid.uuid4())
    assert is_valid_uuid(uuid_str)


def test_is_invalid_uuid():
    assert not is_valid_uuid("blablabla")
