import csv
import uuid
import pytest

from scripts.db.queries import get_analysis_run, get_analysis_run_sample
from scripts.util.metadata import (
    generate_notifications,
    is_valid_uuid,
    inspect_metadata_file,
    link_sample_to_analysis_run,
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
            "good_metadata_illumina_bam.tsv",
            False,
            [],
        ),
        (
            "good_metadata_illumina_bam.tsv",
            True,
            [
                "file_2 for 37a36d1c-5985-4836-87b5-b36bac75d81b not available\n"
                "md5_2 for 37a36d1c-5985-4836-87b5-b36bac75d81b not available",
                "file_2 for 985347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
                "md5_2 for 985347c5-ff6a-454c-ac34-bc353d05dd70 not available",
            ],
        ),
        (
            "bad_metadata.tsv",
            True,
            [
                "sample_id not available",
                'sample_id "#()aadd" is not a UUID',
                "file_1 for 185347c5-ff6a-454c-ac34-bc353d05dd70 not available",
                "file_2 for 27a36d1c-5985-4836-87b5-b36bac75d81b not available",
                "md5_1 for 385347c5-ff6a-454c-ac34-bc353d05dd70 not available",
                "md5_2 for 485347c5-ff6a-454c-ac34-bc353d05dd70 not available",
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

    generate_notifications(passed_samples, passed_path, failed_samples, failed_path)

    for expected_path, expected_samples in data.items():
        assert expected_path.is_file()
        with open(expected_path, "r") as fd:
            loaded_samples = fd.read().splitlines()
        assert expected_samples == loaded_samples


@pytest.mark.parametrize(
    "sample_name, load_missing_samples, expected_exception",
    [
        ("7284954", False, None),
        ("new_sample", True, None),
        (
            "new_sample",
            True,
            ValueError("Sample new_sample not found in the database, but listed in pipeline metadata"),
        ),
    ],
)
def test_link_sample_to_analysis_run(
    db_session, populated_db_session_with_sample, sample_name, load_missing_samples, expected_exception
):

    analysis_run = get_analysis_run(db_session, "just_a_name")

    try:
        link_sample_to_analysis_run(db_session, analysis_run, sample_name, load_missing_samples)

        sample = get_analysis_run_sample(db_session, analysis_run.analysis_run_name, sample_name)
        assert sample is not None
        assert sample.sample_name == sample_name
        assert sample.analysis_run_id == analysis_run.analysis_run_id

    except ValueError as ex:
        assert ex == expected_exception


@pytest.mark.parametrize(
    "metadata_file,samples_with_two_reads,load_missing_samples,valid_samples,invalid_samples",
    [
        (
            "good_metadata_illumina_bam.tsv",
            False,
            True,
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
        ),
        (
            "good_metadata_illumina_bam.tsv",
            False,
            False,
            [],
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
        ),
        (
            "good_metadata_illumina_bam.tsv",
            True,
            False,
            [],
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
        ),
        (
            "bad_metadata.tsv",
            True,
            False,
            [],
            [
                "",
                "#()aadd",
                "185347c5-ff6a-454c-ac34-bc353d05dd70",
                "27a36d1c-5985-4836-87b5-b36bac75d81b",
                "385347c5-ff6a-454c-ac34-bc353d05dd70",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
        ),
    ],
)
def test_inspect_metadata_file(
    db_session,
    populated_db_session_with_sample,
    test_data_path,
    metadata_file,
    samples_with_two_reads,
    load_missing_samples,
    valid_samples,
    invalid_samples,
):
    metadata_path = test_data_path / metadata_file
    analysis_run = get_analysis_run(db_session, "just_a_name")

    samples = inspect_metadata_file(
        db_session, metadata_path, analysis_run, samples_with_two_reads, load_missing_samples
    )

    assert sorted(samples.valid) == sorted(valid_samples)
    assert sorted(samples.invalid) == sorted(invalid_samples)


def test_inspect_metadata_file_exception(
    db_session,
    populated_db_session_with_sample,
    test_data_path,
):
    metadata_path = test_data_path / "good_metadata_illumina_bam.tsv"

    try:
        inspect_metadata_file(db_session, metadata_path, None, False, False)
    except ValueError as err:
        assert str(err) == "Cannot retrieve samples for a non existing analysis run."
