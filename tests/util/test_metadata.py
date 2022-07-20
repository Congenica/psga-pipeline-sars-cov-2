import csv
import uuid
import pytest

from scripts.util.metadata import (
    generate_notifications,
    is_valid_uuid,
    inspect_metadata_file,
    validate_and_normalise_row,
)


def test_is_valid_uuid():
    uuid_str = str(uuid.uuid4())
    assert is_valid_uuid(uuid_str)


def test_is_invalid_uuid():
    assert not is_valid_uuid("blablabla")


@pytest.mark.parametrize(
    "metadata_file,samples_with_two_reads,expected_exceptions",
    [
        (
            "good_metadata_illumina_bam.csv",
            False,
            [],
        ),
        (
            "good_metadata_illumina_bam.csv",
            True,
            [
                "file_2 for 37a36d1c-5985-4836-87b5-b36bac75d81b not available",
                "file_2 for 985347c5-ff6a-454c-ac34-bc353d05dd70 not available",
            ],
        ),
        (
            "bad_metadata.csv",
            True,
            [
                "sample_id not available",
                'sample_id "#()aadd" is not a UUID',
                "file_1 for 185347c5-ff6a-454c-ac34-bc353d05dd70 not available",
                "file_2 for 27a36d1c-5985-4836-87b5-b36bac75d81b not available",
            ],
        ),
    ],
)
def test_validate_and_normalise_row(
    test_data_path,
    metadata_file,
    samples_with_two_reads,
    expected_exceptions,
):
    metadata_path = test_data_path / metadata_file
    exceptions = []
    with open(metadata_path, "r") as metadata_fd:
        reader = csv.DictReader(metadata_fd, delimiter=",")
        for row in reader:
            try:
                row = validate_and_normalise_row(row, samples_with_two_reads)
            except ValueError as e:
                exceptions.append(str(e))
    assert sorted(exceptions) == sorted(expected_exceptions)


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


@pytest.mark.parametrize(
    "metadata_file,samples_with_two_reads,valid_samples,invalid_samples",
    [
        (
            "good_metadata_illumina_bam.csv",
            False,
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
        ),
        (
            "bad_metadata.csv",
            True,
            ["385347c5-ff6a-454c-ac34-bc353d05dd70"],
            [
                "",
                "#()aadd",
                "185347c5-ff6a-454c-ac34-bc353d05dd70",
                "27a36d1c-5985-4836-87b5-b36bac75d81b",
            ],
        ),
    ],
)
def test_inspect_metadata_file(
    test_data_path,
    metadata_file,
    samples_with_two_reads,
    valid_samples,
    invalid_samples,
):
    metadata_path = test_data_path / metadata_file

    samples = inspect_metadata_file(metadata_path, samples_with_two_reads)

    assert sorted(samples.valid) == sorted(valid_samples)
    assert sorted(samples.invalid) == sorted(invalid_samples)
