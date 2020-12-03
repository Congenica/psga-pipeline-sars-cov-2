import csv
from pathlib import Path
from typing import List, Tuple, Any

import click
import requests
from requests.exceptions import ConnectTimeout, MissingSchema
from sqlalchemy import func

from db.database import session_handler
from db.models import Sample

STRAIN_LEVEL_AND_GLOBAL_CONTEXT_REPORT_HEADERS = ("Lineage", "Description", "Count", "Percent", "Global Count")


def fetch_csv(url: str, delimiter: str, fallback_dir: str):
    """
    Function to try fetching latest data, if fetch fails, looks up local files specified in fallback_dir
    :param url: url to try fetching the latest csv file
    :param delimiter: one-character string to use as the field separator
    :param fallback_dir: fallback if latest data fetch fails, will look for latest local data by filename from url
    :return: Tuple[List[Any, ...]]
    for notes the return is converted into dictionary, but mypy doesn't understand the structure,
    so return typing is omitted
    """

    def parse_csv(content):
        content = csv.reader(content, delimiter=delimiter)
        return tuple(content)

    try:
        with requests.Session() as s:
            response = s.get(url)
            if not response.ok:
                raise ConnectionError()
            decoded_content = response.content.decode("utf-8")
            return parse_csv(decoded_content.splitlines())
    except (ConnectionError, ConnectTimeout, MissingSchema):
        filename = url.split("/")[-1]
        path = list(Path(fallback_dir).rglob(f"{filename}"))[0]
        with open(path, "r") as csv_file:
            return parse_csv(csv_file)


def get_strain_level_and_global_context_report_data(
    lineage_notes_url: str, metadata_url: str, pangolearn_dir: str
) -> List[Tuple[Any, ...]]:
    """
    Method to fetch report data from local database and external data sources.
    For more details on param usage lookup fetch_csv docstring
    :param str lineage_notes_url: url to latest pangoLEARN lineage notes
    :param str metadata_url: url to latest pangoLEARN metadata
    :param str pangolearn_dir: fallback pangoLEARN directory with local lineage notes and metadata
    :return List[Tuple[Any, ...]]: lines for the csv report containing columns:
    Lineage, Description (from lineage notes), count, percent, global count (from metadata)
    order corresponds with STRAIN_LEVEL_AND_GLOBAL_CONTEXT_REPORT_HEADERS
    """
    # fetch database data
    with session_handler() as session:
        lineage_count = session.query(Sample.pangolin_lineage, func.count("*")).group_by(Sample.pangolin_lineage).all()
    # fetch lineage notes
    notes = fetch_csv(lineage_notes_url, "\t", pangolearn_dir)
    # fetch metadata
    meta_data = fetch_csv(metadata_url, ",", pangolearn_dir)
    # adapt report values
    result: List[Tuple[Any, ...]] = [STRAIN_LEVEL_AND_GLOBAL_CONTEXT_REPORT_HEADERS]
    notes = dict(notes)
    total = sum(x[1] for x in lineage_count)
    global_total = 0
    for lineage in lineage_count:
        if lineage[0] is None:
            continue
        global_count = sum(1 for x in meta_data if x[-1] == lineage[0])
        global_total += global_count
        percentage = "{:.0%}".format(lineage[1] / total)
        result.append(
            (
                lineage[0],
                notes[lineage[0]],
                lineage[1],
                percentage,
                global_count,
            )
        )
    # footer
    result.append(
        (
            "",
            "",
            total,
            "",
            global_total,
        )
    )
    return result


reports = {
    "strain_level_and_global_context": get_strain_level_and_global_context_report_data,
}


@click.command()
@click.option(
    "--output",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="output report file path, where the report will be written to",
)
@click.option(
    "--report",
    type=click.Choice(["strain_level_and_global_context"]),
    required=True,
    help="",
)
@click.option(
    "--lineage-notes-url",
    type=click.STRING,
)
@click.option(
    "--metadata-url",
    type=click.STRING,
)
@click.option(
    "--pangolearn-dir",
    type=click.Path(dir_okay=True, readable=True),
)
def generate_report(output: str, report: str, **kwargs) -> None:
    """
    Generate a csv report from the database and external sources
    """
    samples = reports[report](**kwargs)

    with open(output, "w") as out_file:
        csv_writer = csv.writer(out_file)
        csv_writer.writerows(samples)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_report()
