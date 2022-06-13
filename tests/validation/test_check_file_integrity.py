import pytest
from click.testing import CliRunner

from scripts.validation.check_file_integrity import check_file_integrity, FileIntegrityError


@pytest.mark.parametrize(
    "input_path,expected_md5,exit_code,exception",
    [
        (
            "good_metadata_illumina_fastq.tsv",
            "12e5fd2bc82a787d40c4e014e67a0295",
            0,
            None,
        ),
        (
            "good_metadata_illumina_fastq.tsv",
            "fake_md5",
            1,
            FileIntegrityError(
                "Integrity check for file: good_metadata_illumina_fastq.tsv FAILED. "
                + "Expected: fake_md5, computed: 12e5fd2bc82a787d40c4e014e67a0295"
            ),
        ),
    ],
)
def test_check_metadata(
    test_data_path,
    input_path,
    expected_md5,
    exit_code,
    exception,
):

    cmd_config = [
        "--analysis-run-name",
        "just_a_run",
        "--sample-name",
        "what_a_sample",
        "--input-path",
        test_data_path / input_path,
        "--expected-md5",
        expected_md5,
    ]

    rv = CliRunner().invoke(
        check_file_integrity,
        cmd_config,
    )

    assert rv.exit_code == exit_code

    if exit_code == 1:
        assert rv.exception == exception
