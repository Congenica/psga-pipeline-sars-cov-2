import uuid

from scripts.util.metadata import (
    generate_notifications,
    is_valid_uuid,
)


def test_is_valid_uuid():
    uuid_str = str(uuid.uuid4())
    assert is_valid_uuid(uuid_str)


def test_is_invalid_uuid():
    assert not is_valid_uuid("blablabla")


def test_generate_notifications(tmp_path):
    analysis_run = "analysis_run1"
    passed_path = tmp_path / "passed"
    passed_samples = ["v1", "v2", "v3"]
    failed_path = tmp_path / "failed"
    failed_samples = ["x1", "x2"]
    data = {
        passed_path: passed_samples,
        failed_path: failed_samples,
    }

    for f in data:
        assert not f.is_file()

    generate_notifications(analysis_run, passed_samples, passed_path, failed_samples, failed_path)

    for expected_path, expected_samples in data.items():
        assert expected_path.is_file()
        with open(expected_path, "r") as fd:
            loaded_samples = fd.read().splitlines()
        assert expected_samples == loaded_samples
