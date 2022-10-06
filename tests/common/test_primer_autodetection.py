import pytest

from click.testing import CliRunner

from scripts.common.primer_autodetection import (
    select_primer,
    write_primer_data,
    write_selected_primer,
    generate_primer_autodetection_output_files,
    primer_autodetection,
    PRIMER_AUTODETECTION_PRIMER_COL,
    PRIMER_AUTODETECTION_PRIMER_SCORE_COL,
    PRIMER_AUTODETECTION_SAMPLE_ID_COL,
    PRIMER_DETECTION_SUFFIX,
    PRIMER_DATA_SUFFIX,
    UNKNOWN,
)
from tests.utils_tests import assert_csvs_are_equal, assert_files_are_equal

DEFAULT_PRIMER = "ARTIC_V4-1"


def assert_primer_detection(sample_id, tmp_path, input_path):
    output_file = f"{sample_id}{PRIMER_DETECTION_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_csvs_are_equal(output_path, expected_output_path, PRIMER_AUTODETECTION_PRIMER_COL)


def assert_primer_data(sample_id, tmp_path, input_path):
    output_file = f"{sample_id}{PRIMER_DATA_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_csvs_are_equal(output_path, expected_output_path, PRIMER_AUTODETECTION_SAMPLE_ID_COL)


def assert_selected_primer_file(sample_id, test_data_path, tmp_path, found_dir):
    primer_file = f"{sample_id}_primer.txt"
    expected_output_path = test_data_path / "primer_autodetection" / found_dir / primer_file
    output_path = tmp_path / primer_file
    assert_files_are_equal(output_path, expected_output_path)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input,expected_score,expected_detected_primer,expected_selected_primer",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "ARTIC_V3", 423, "ARTIC_V4-1", "ARTIC_V3"),
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", UNKNOWN, 423, "ARTIC_V4-1", "ARTIC_V4-1"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "ARTIC_V3", 0, "none", "ARTIC_V3"),
        # test the hack
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", UNKNOWN, 0, "none", DEFAULT_PRIMER),
    ],
)
def test_select_primer(
    tmp_path,
    test_data_path,
    found_dir,
    sample_id,
    primer_input,
    expected_score,
    expected_detected_primer,
    expected_selected_primer,
):
    input_path = test_data_path / "primer_autodetection" / found_dir
    detected_primer_df, selected_primer = select_primer(input_path, tmp_path, sample_id, primer_input)

    assert expected_detected_primer == detected_primer_df[PRIMER_AUTODETECTION_PRIMER_COL]
    assert expected_score == detected_primer_df[PRIMER_AUTODETECTION_PRIMER_SCORE_COL]
    assert expected_selected_primer == selected_primer
    assert_primer_detection(sample_id, tmp_path, input_path)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown"),
    ],
)
def test_write_primer_data(tmp_path, test_data_path, found_dir, sample_id, primer_input):
    input_path = test_data_path / "primer_autodetection" / found_dir
    detected_primer_df, _ = select_primer(input_path, tmp_path, sample_id, primer_input)
    write_primer_data(tmp_path, detected_primer_df, sample_id, primer_input)
    assert_primer_data(sample_id, tmp_path, input_path)


@pytest.mark.parametrize(
    "found_dir,sample_id,detected_primer",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "ARTIC_V4-1"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", DEFAULT_PRIMER),
    ],
)
def test_write_selected_primer(tmp_path, test_data_path, found_dir, sample_id, detected_primer):
    write_selected_primer(tmp_path, sample_id, detected_primer)
    assert_selected_primer_file(sample_id, test_data_path, tmp_path, found_dir)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown"),
    ],
)
def test_generate_primer_autodetection_output_files(tmp_path, test_data_path, found_dir, sample_id, primer_input):
    input_path = test_data_path / "primer_autodetection" / found_dir
    generate_primer_autodetection_output_files(input_path, tmp_path, sample_id, primer_input)
    assert_primer_detection(sample_id, tmp_path, input_path)
    assert_primer_data(sample_id, tmp_path, input_path)
    assert_selected_primer_file(sample_id, test_data_path, tmp_path, found_dir)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown"),
    ],
)
def test_primer_autodetection(tmp_path, test_data_path, found_dir, sample_id, primer_input):
    input_path = test_data_path / "primer_autodetection" / found_dir

    rv = CliRunner().invoke(
        primer_autodetection,
        [
            "--input-path",
            input_path,
            "--output-path",
            tmp_path,
            "--sample-id",
            sample_id,
            "--primer-input",
            primer_input,
        ],
    )
    assert rv.exit_code == 0

    assert_primer_detection(sample_id, tmp_path, input_path)
    assert_primer_data(sample_id, tmp_path, input_path)
    assert_selected_primer_file(sample_id, test_data_path, tmp_path, found_dir)
