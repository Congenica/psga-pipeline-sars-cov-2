import pytest
from click.testing import CliRunner

from scripts.common.format_genotyping import (
    process_typing,
    format_genotyping,
    TYPING_SAMPLE_ID_COL,
)
from tests.utils_tests import assert_dataframes_are_equal


@pytest.mark.parametrize(
    "input_file,input_sample_id_col,input_type_col,expected_typing_output_csv",
    [
        (
            "one_sample_two_types.csv",
            "sampleID",
            "type",
            "expected_output_one_sample_two_types.csv",
        ),
        (
            "two_samples_one_type.csv",
            "sampleID",
            "type",
            "expected_output_two_samples_one_type.csv",
        ),
    ],
)
def test_process_typing(
    tmp_path, test_data_path, input_file, input_sample_id_col, input_type_col, expected_typing_output_csv
):
    input_path = test_data_path / "typing" / input_file
    expected_typing_output_path = test_data_path / "typing" / expected_typing_output_csv
    output_path = tmp_path / "typing_output.csv"
    process_typing(input_path, output_path, input_sample_id_col, input_type_col)
    assert_dataframes_are_equal(output_path, expected_typing_output_path, TYPING_SAMPLE_ID_COL)


@pytest.mark.parametrize(
    "input_file,input_sample_id_col,input_type_col,expected_typing_output_csv",
    [
        (
            "one_sample_two_types.csv",
            "sampleID",
            "type",
            "expected_output_one_sample_two_types.csv",
        ),
        (
            "two_samples_one_type.csv",
            "sampleID",
            "type",
            "expected_output_two_samples_one_type.csv",
        ),
    ],
)
def test_typing(tmp_path, test_data_path, input_file, input_sample_id_col, input_type_col, expected_typing_output_csv):
    input_path = test_data_path / "typing" / input_file
    expected_typing_output_path = test_data_path / "typing" / expected_typing_output_csv
    output_path = tmp_path / "typing_output.csv"
    rv = CliRunner().invoke(
        format_genotyping,
        [
            "--input-path",
            input_path,
            "--output-csv-path",
            output_path,
            "--input-sample-id-col",
            input_sample_id_col,
            "--input-type-col",
            input_type_col,
        ],
    )
    assert rv.exit_code == 0
    assert_dataframes_are_equal(output_path, expected_typing_output_path, TYPING_SAMPLE_ID_COL)
