from pathlib import Path
import pandas as pd


def _refine_df(config: dict, df: pd.DataFrame) -> pd.DataFrame:
    """
    Call common functions for refining the dataframe
    """

    for col in ["sample_name_column", "columns_to_validate", "columns_to_round"]:
        if col not in config:
            raise ValueError(f"Error in the validation config. {col} not specified in config")

    sample_name_column = config["sample_name_column"]
    columns_to_validate = config["columns_to_validate"]

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


def load_data_from_csv(config: dict, csv_path: Path) -> pd.DataFrame:
    """
    Load the CSV content to a Pandas dataframe, performing basic validation
    """
    df = pd.read_csv(csv_path)

    return _refine_df(config, df)


def get_file_paths(root: Path) -> list[str]:
    """
    Return a list of file paths in root. The search is recursive.
    """
    file_list = []
    for x in root.iterdir():
        if x.is_file():
            file_list.append(str(x))
        elif x.is_dir():
            file_list.extend(get_file_paths(x))
    return file_list
