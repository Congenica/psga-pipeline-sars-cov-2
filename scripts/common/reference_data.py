from csv import DictReader
from itertools import chain
import pathlib
from typing import Iterator
import click


class ReferenceDataNotFoundError(Exception):
    """
    Exception raised when the given Reference Data set could not
    be found in the Reference Data CSV file.
    """


def normalise_csv(file_lines: Iterator[str]) -> Iterator[str]:
    """
    Returns the same iterator with the CSV file lines but after
    converting the first row (headers) to lowercase.

    At the time of implementing, the headers of the Reference Data
    CSV file looked like "NAME", and "LOCATION". This normalisation
    makes the script more resilient.
    """
    return chain([next(file_lines).lower()], file_lines)


def reference_data_csv_to_dict(reference_data_file: pathlib.Path) -> None:
    """
    :param reference_data_file: The path to the Reference Data CSV file.
    :returns: Iterates over the Reference Data CSV file and returns a dictionary where the
        keys are the names of the Reference Data sets, and the values are the locations.
    """
    reference_data_sets: dict[str, pathlib.Path] = {}
    with reference_data_file.open() as fd:
        reader = DictReader(normalise_csv(file_lines=fd))
        for row in reader:
            name = row["name"]
            location = pathlib.Path(row["location"])
            reference_data_sets[name] = location
    return reference_data_sets


@click.group()
def cli():
    pass


@click.command()
@click.argument(
    "reference_data_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
@click.argument(
    "name",
    type=str,
)
def get_location(reference_data_file: pathlib.Path, name: str) -> None:
    """
    Using the Reference Data CSV file, this command returns the location of
    the Reference Data set identified by the given name.

    REFERENCE_DATA_FILE is the path to the Reference Data CSV file.

    NAME is the string identifying a Reference Data set.
    """
    reference_data_sets = reference_data_csv_to_dict(reference_data_file=reference_data_file)
    if name not in reference_data_sets:
        raise ReferenceDataNotFoundError(
            f'The Reference Data with name "{name}" was not found '
            f'in the CSV file "{reference_data_file}".'
        )
    location: pathlib.Path = reference_data_sets[name]
    print(location.absolute())


cli.add_command(get_location)


if __name__ == "__main__":
    cli()
