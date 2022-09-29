import os
import pytest

import pandas as pd

from scripts.util.metadata import SAMPLE_ID
from jenkins.loading import get_file_paths, load_data_from_csv
from tests.jenkins.util import create_paths


@pytest.mark.parametrize(
    "csv_file,config,exc",
    [
        (
            "merged_output.csv",
            {},
            "sample_name_column not specified in config",
        ),
        (
            "merged_output.csv",
            {
                "sample_name_column": SAMPLE_ID,
            },
            "columns_to_validate not specified in config",
        ),
        (
            "merged_output.csv",
            {
                "sample_name_column": SAMPLE_ID,
                "columns_to_validate": [SAMPLE_ID],
            },
            "columns_to_round not specified in config",
        ),
        (
            "merged_output.csv",
            {
                "sample_name_column": SAMPLE_ID,
                "columns_to_validate": [SAMPLE_ID],
                "columns_to_round": [],
            },
            None,
        ),
        (
            "merged_output.csv",
            {
                "sample_name_column": "fake_sample_id",
                "columns_to_validate": [SAMPLE_ID],
                "columns_to_round": [],
            },
            "Make sure that fake_sample_id is in 'columns_to_validate'",
        ),
        (
            "merged_output.csv",
            {
                "sample_name_column": SAMPLE_ID,
                "columns_to_validate": ["fake_sample_id"],
                "columns_to_round": [],
            },
            f"Make sure that {SAMPLE_ID} is in 'columns_to_validate'",
        ),
    ],
)
def test_load_data_from_csv(tmp_path, test_data_path, csv_file, config, exc):
    csv_path = test_data_path / "jenkins" / "validation" / csv_file
    # simplify this file for testing
    df = pd.read_csv(csv_path)
    # with the following op, the dataframe becomes a series
    df = df[SAMPLE_ID]

    df_path = tmp_path / csv_file
    df.to_csv(df_path)

    if exc:
        with pytest.raises(ValueError, match=exc):
            load_data_from_csv(config, csv_path)
    else:
        df_loaded = load_data_from_csv(config, csv_path)
        assert set(df.tolist()) == set(df_loaded[SAMPLE_ID].tolist())


@pytest.mark.parametrize(
    "path_list",
    [
        [],
        ["a/b/c1", "a/b/c2", "d/e", "f"],
        ["a", "b", "c"],
    ],
)
def test_get_file_paths(tmp_path, path_list):
    paths = [os.path.join(tmp_path, f) for f in path_list]
    create_paths(paths)
    fetched_paths = get_file_paths(tmp_path)
    assert set(fetched_paths) == set(paths)
