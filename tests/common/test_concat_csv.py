from pathlib import Path
import pytest
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from app.scripts.concat_csv import concat_csv, concat
from app.scripts.util.metadata import SAMPLE_ID
from tests.utils_tests import assert_csvs_are_equal


@pytest.mark.jira(identifier="841d3e5c-65ea-4db2-8247-713dc550ac32", confirms="PSG-3621")
def test_concat_no_file_found_error(tmp_path: Path):
    output_csv_path = Path(tmp_path / "concat.csv")
    with pytest.raises(FileNotFoundError, match="No matching file was found"):
        concat(input_path=tmp_path, output_csv_path=output_csv_path)


@pytest.mark.jira(identifier="edd1fee1-2b71-4c4f-bf45-c8de604610cc", confirms="PSG-3621")
def test_concat_column_not_found_error(tmp_path: Path, concat_csv_data_path: Path):
    output_csv_path = Path(tmp_path / "concat.csv")
    with pytest.raises(ValueError, match="Column 'fake' needed for sorting"):
        concat(input_path=concat_csv_data_path, output_csv_path=output_csv_path, sortby_col="fake")


@pytest.mark.jira(identifier="e9b930ec-9b72-401c-8933-bc17d2fd5238", confirms="PSG-3621")
def test_concat(tmp_path: Path, concat_csv_data_path: Path):
    output_csv_path = Path(tmp_path / "concat.csv")
    expected_output_csv_path = Path(concat_csv_data_path / "expected_output" / "result.csv")

    assert not output_csv_path.is_file()
    concat_df = concat(input_path=concat_csv_data_path, output_csv_path=output_csv_path, sortby_col=SAMPLE_ID)

    assert output_csv_path.is_file()

    assert_csvs_are_equal(output_csv_path, expected_output_csv_path, SAMPLE_ID)

    df_exp = pd.read_csv(expected_output_csv_path)

    assert_frame_equal(
        concat_df.reset_index(drop=True),
        df_exp.reset_index(drop=True),
    )


@pytest.mark.jira(identifier="7aaa925b-e92a-49a7-b0df-0a781b25d9c5", confirms="PSG-3621")
def test_concat_csv(tmp_path: Path, concat_csv_data_path: Path):
    output_csv_path = Path(tmp_path / "concat.csv")
    expected_output_csv_path = Path(concat_csv_data_path / "expected_output" / "result.csv")

    rv = CliRunner().invoke(
        concat_csv,
        [
            "--input-path",
            concat_csv_data_path,
            "--output-csv-path",
            output_csv_path,
            "--sortby-col",
            SAMPLE_ID,
        ],
    )

    assert rv.exit_code == 0

    assert_csvs_are_equal(output_csv_path, expected_output_csv_path, SAMPLE_ID)
