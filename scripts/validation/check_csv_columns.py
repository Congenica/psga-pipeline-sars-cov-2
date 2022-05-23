from typing import Set
from click import ClickException


def check_csv_columns(reader_fieldnames: Set[str], expected_fieldnames: Set[str]) -> None:
    """
    Check whether the expected fieldnames are present in the reader fieldnames
    """
    if not expected_fieldnames.issubset(reader_fieldnames):
        err = (
            "Unexpected headers, got:\n"
            + ", ".join(reader_fieldnames)
            + "\n, but expect at least \n"
            + ", ".join(expected_fieldnames)
        )
        raise ClickException(err)
