from pathlib import Path
from typing import Dict
import json
import pandas as pd
from pandas.testing import assert_frame_equal


def assert_files_are_equal(file1: Path, file2: Path) -> None:
    with open(file1, "r") as f1:
        with open(file2, "r") as f2:
            content_1 = f1.read()
            content_2 = f2.read()
            assert content_1 == content_2


def assert_dataframes_are_equal(file1: Path, file2: Path, sortby_col: str) -> None:
    df = pd.read_csv(file1)
    df_exp = pd.read_csv(file2)

    # sort rows by sortby_col
    df_sorted = df.sort_values(by=[sortby_col])
    df_exp_sorted = df_exp.sort_values(by=[sortby_col])

    # sort columns by column name
    df_sorted = df_sorted.reindex(columns=sorted(df_sorted.columns))
    df_exp_sorted = df_exp_sorted.reindex(columns=sorted(df_exp_sorted.columns))
    assert_frame_equal(df_sorted, df_exp_sorted)


def json_to_dict(input_path: Path) -> Dict:
    with open(input_path, "r", encoding="utf-8") as jsonfile:
        return json.load(jsonfile)


def dict_to_string(input_dict: Dict) -> str:
    # sort all keys, including nested dictionaries and lists
    return json.dumps(input_dict, sort_keys=True)


def json_to_string(input_path: Path) -> str:
    return dict_to_string(json_to_dict(input_path))


def assert_jsons_are_equal(file1: Path, file2: Path) -> None:
    assert json_to_string(file1) == json_to_string(file2)
