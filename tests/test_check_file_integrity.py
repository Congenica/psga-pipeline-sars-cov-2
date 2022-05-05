import pytest
from click.testing import CliRunner

from scripts.check_file_integrity import check_file_integrity, FileIntegrityError


@pytest.mark.parametrize(
    "input_path,expected_md5,exit_code,exception",
    [
        (
            "good_metadata.tsv",
            "91ec15a673c2692949847111b0366bc3",
            0,
            None,
        ),
        (
            "good_metadata.tsv",
            "fake_md5",
            1,
            FileIntegrityError(
                "Integrity check for file: good_metadata.tsv FAILED. "
                + "Expected: fake_md5, computed: 91ec15a673c2692949847111b0366bc3"
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
