from click import ClickException


def check_csv_columns(reader_fieldnames: set[str], expected_fieldnames: set[str]) -> None:
    """
    Check whether the expected fieldnames are present in the reader fieldnames
    """
    if expected_fieldnames != reader_fieldnames:
        err = (
            "Unexpected headers, got:\n"
            + ", ".join(sorted(reader_fieldnames))
            + "\n, but expect at least \n"
            + ", ".join(sorted(expected_fieldnames))
        )
        raise ClickException(err)
