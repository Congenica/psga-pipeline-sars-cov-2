from pathlib import Path
import pytest
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.common.concat_csv import concat_csv, concat
from scripts.util.metadata import SAMPLE_ID
from tests.utils_tests import assert_csvs_are_equal


def test_concat_no_file_found_error(tmp_path):
    output_csv_path = Path(tmp_path / "concat.csv")
    with pytest.raises(FileNotFoundError, match="No matching file was found"):
        concat(input_path=tmp_path, output_csv_path=output_csv_path)


def test_concat_column_not_found_error(tmp_path, test_data_path):
    input_path = Path(test_data_path / "concat_csv")
    output_csv_path = Path(tmp_path / "concat.csv")
    with pytest.raises(ValueError, match="Column 'fake' needed for sorting"):
        concat(input_path=input_path, output_csv_path=output_csv_path, sortby_col="fake")


def test_concat(tmp_path, test_data_path):
    input_path = Path(test_data_path / "concat_csv")
    output_csv_path = Path(tmp_path / "concat.csv")
    expected_output_csv_path = Path(test_data_path / "concat_csv" / "expected_output" / "result.csv")

    assert not output_csv_path.is_file()
    concat_df = concat(input_path=input_path, output_csv_path=output_csv_path, sortby_col=SAMPLE_ID)

    assert output_csv_path.is_file()

    assert_csvs_are_equal(output_csv_path, expected_output_csv_path, SAMPLE_ID)

    df_exp = pd.read_csv(expected_output_csv_path)

    assert_frame_equal(
        concat_df.reset_index(drop=True),
        df_exp.reset_index(drop=True),
    )


def test_concat_csv(tmp_path, test_data_path):

    input_path = Path(test_data_path / "concat_csv")
    output_csv_path = Path(tmp_path / "concat.csv")
    expected_output_csv_path = Path(test_data_path / "concat_csv" / "expected_output" / "result.csv")

    rv = CliRunner().invoke(
        concat_csv,
        [
            "--input-path",
            input_path,
            "--output-csv-path",
            output_csv_path,
            "--sortby-col",
            SAMPLE_ID,
        ],
    )

    assert rv.exit_code == 0

    assert_csvs_are_equal(output_csv_path, expected_output_csv_path, SAMPLE_ID)
