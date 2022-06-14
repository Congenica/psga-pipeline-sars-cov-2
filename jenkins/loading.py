from pathlib import Path
from typing import Dict
import pandas as pd


def refine_df(config: Dict, df: pd.DataFrame) -> pd.DataFrame:
    """
    Call common functions for refining the dataframe
    """
    sample_name_column = config["sample_name_column"]
    columns_to_validate = config["columns_to_validate"].values()
    if sample_name_column not in columns_to_validate:
        raise ValueError(
            f"Error in the validation config. Make sure that {sample_name_column} is in 'columns_to_validate'"
        )
    if sample_name_column not in df.columns:
        raise ValueError(
            f"The column {sample_name_column} was not found in the data frame. Impossible to retrieve samples"
        )

    # select the columns to validate and sort them by sample_name_column
    df = df[columns_to_validate].sort_values(by=[sample_name_column])

    return df


def load_data_from_csv(config: Dict, csv_path: Path) -> pd.DataFrame:
    """
    Load the CSV content to a Pandas dataframe, performing basic validation
    """
    df = pd.read_csv(csv_path)

    # rename columns so that they match the names used in the DB
    columns_mapping = config["columns_to_validate"]
    df = df.rename(columns=columns_mapping)

    return refine_df(config, df)
