from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal


def assert_files_are_equal(file1: Path, file2: Path):
    with open(file1, "r") as f1:
        with open(file2, "r") as f2:
            content_1 = f1.read()
            content_2 = f2.read()
            assert content_1 == content_2


def assert_dataframes_are_equal(file1: Path, file2: Path, sortby_col: str):
    df = pd.read_csv(file1)
    df_exp = pd.read_csv(file2)

    df_sorted = df.sort_values(by=[sortby_col])
    df_exp_sorted = df_exp.sort_values(by=[sortby_col])
    assert_frame_equal(
        df_sorted.reset_index(drop=True),
        df_exp_sorted.reset_index(drop=True),
    )
