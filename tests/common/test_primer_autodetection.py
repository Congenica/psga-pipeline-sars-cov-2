import pytest

from click.testing import CliRunner

from scripts.common.primer_autodetection import (
    select_primer,
    calculate_primer_qc,
    write_primer_data,
    write_selected_primer,
    generate_primer_autodetection_output_files,
    primer_autodetection,
    PRIMER_AUTODETECTION_PRIMER_COL,
    PRIMER_AUTODETECTION_PRIMER_SCORE_COL,
    PRIMER_AUTODETECTION_SAMPLE_ID_COL,
    PRIMER_QC_PASS,
    PRIMER_QC_FAIL,
    PRIMER_DETECTION_SUFFIX,
    PRIMER_DATA_SUFFIX,
    RNAME_COL,
)
from tests.utils_tests import assert_dataframes_are_equal, assert_files_are_equal


def assert_primer_detection(sample_id, tmp_path, input_path):
    output_file = f"{sample_id}{PRIMER_DETECTION_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_dataframes_are_equal(output_path, expected_output_path, RNAME_COL)


def assert_primer_data(sample_id, tmp_path, input_path):
    output_file = f"{sample_id}{PRIMER_DATA_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_dataframes_are_equal(output_path, expected_output_path, PRIMER_AUTODETECTION_SAMPLE_ID_COL)


def assert_selected_primer_file(sample_id, primer_qc, test_data_path, tmp_path, qc_dir):
    primer_file = f"{sample_id}_primer_{primer_qc}.txt"
    expected_output_path = test_data_path / "primer_autodetection" / qc_dir / primer_file
    output_path = tmp_path / primer_file
    assert_files_are_equal(output_path, expected_output_path)


@pytest.mark.parametrize(
    "qc_dir,sample_id,primer,score",
    [
        ("pass", "9729bce7-f0a9-4617-b6e0-6145307741d1", "ARTIC_V4-1", 4),
        ("fail", "a0446f6f-7d24-478c-8d92-7c77036930d8", "none", 0),
    ],
)
def test_select_primer(tmp_path, test_data_path, qc_dir, sample_id, primer, score):
    input_path = test_data_path / "primer_autodetection" / qc_dir
    detected_primer_df = select_primer(input_path, tmp_path, sample_id)
    assert primer == detected_primer_df[PRIMER_AUTODETECTION_PRIMER_COL]
    assert score == detected_primer_df[PRIMER_AUTODETECTION_PRIMER_SCORE_COL]
    assert_primer_detection(sample_id, tmp_path, input_path)


@pytest.mark.parametrize(
    "primer_detected, primer_score",
    [
        ("none", 4),
        ("ARTIC_V4", 0),
    ],
)
def test_calculate_primer_qc_error(primer_detected, primer_score):
    with pytest.raises(ValueError, match="Error in computing primer QC"):
        calculate_primer_qc("fake", primer_detected, primer_score)


@pytest.mark.parametrize(
    "primer_input, primer_detected, primer_score, expected_qc",
    [
        # user expects sample to have primer sequences but does not know the primer
        ("unknown", "ARTIC_V4", 4, PRIMER_QC_PASS),
        ("unknown", "none", 0, PRIMER_QC_FAIL),
        # user knows that sample has no primer sequences
        ("none", "ARTIC_V4", 4, PRIMER_QC_FAIL),
        ("none", "none", 0, PRIMER_QC_PASS),
        # user thinks to know the primer
        ("ARTIC_V4", "ARTIC_V4", 4, PRIMER_QC_PASS),
        ("ARTIC_V3", "ARTIC_V4", 4, PRIMER_QC_FAIL),
    ],
)
def test_calculate_primer_qc(primer_input, primer_detected, primer_score, expected_qc):
    assert calculate_primer_qc(primer_input, primer_detected, primer_score) == expected_qc


@pytest.mark.parametrize(
    "qc_dir,sample_id,primer_input,primer_qc",
    [
        ("pass", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown", PRIMER_QC_PASS),
        ("fail", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown", PRIMER_QC_FAIL),
    ],
)
def test_write_primer_data(tmp_path, test_data_path, qc_dir, sample_id, primer_input, primer_qc):
    input_path = test_data_path / "primer_autodetection" / qc_dir
    detected_primer_df = select_primer(input_path, tmp_path, sample_id)
    write_primer_data(tmp_path, detected_primer_df, sample_id, primer_qc, primer_input)
    assert_primer_data(sample_id, tmp_path, input_path)


@pytest.mark.parametrize(
    "qc_dir,sample_id,primer_input,detected_primer,primer_qc",
    [
        ("pass", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown", "ARTIC_V4-1", PRIMER_QC_PASS),
        ("fail", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown", "none", PRIMER_QC_FAIL),
    ],
)
def test_write_selected_primer(tmp_path, test_data_path, qc_dir, sample_id, primer_input, detected_primer, primer_qc):
    write_selected_primer(tmp_path, sample_id, primer_input, detected_primer, primer_qc)
    assert_selected_primer_file(sample_id, primer_qc, test_data_path, tmp_path, qc_dir)


@pytest.mark.parametrize(
    "qc_dir,sample_id,primer_input,primer_qc",
    [
        ("pass", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown", PRIMER_QC_PASS),
        ("fail", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown", PRIMER_QC_FAIL),
    ],
)
def test_generate_primer_autodetection_output_files(
    tmp_path, test_data_path, qc_dir, sample_id, primer_input, primer_qc
):
    input_path = test_data_path / "primer_autodetection" / qc_dir
    generate_primer_autodetection_output_files(input_path, tmp_path, sample_id, primer_input)
    assert_primer_detection(sample_id, tmp_path, input_path)
    assert_primer_data(sample_id, tmp_path, input_path)
    assert_selected_primer_file(sample_id, primer_qc, test_data_path, tmp_path, qc_dir)


@pytest.mark.parametrize(
    "qc_dir,sample_id,primer_input,primer_qc",
    [
        ("pass", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown", PRIMER_QC_PASS),
        ("fail", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown", PRIMER_QC_FAIL),
    ],
)
def test_primer_autodetection(tmp_path, test_data_path, qc_dir, sample_id, primer_input, primer_qc):
    input_path = test_data_path / "primer_autodetection" / qc_dir

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
    assert_selected_primer_file(sample_id, primer_qc, test_data_path, tmp_path, qc_dir)
