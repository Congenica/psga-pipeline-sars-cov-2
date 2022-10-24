import filecmp
from pathlib import Path
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.util.concat_csv import concat_csv
from scripts.util.metadata import SAMPLE_ID


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

    # check data frames are identical, catching any difference
    df = pd.read_csv(output_csv_path)
    df_exp = pd.read_csv(expected_output_csv_path)

    calculated_cols = df.columns
    expected_cols = df_exp.columns
    assert sorted(calculated_cols) == sorted(expected_cols)

    assert_frame_equal(
        df.reset_index(drop=True),
        df_exp.reset_index(drop=True),
    )

    # check actual csv files are identical
    assert filecmp.cmp(output_csv_path, expected_output_csv_path)
