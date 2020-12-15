import csv
from pathlib import Path
from typing import List, Tuple, Any, Dict, Callable, DefaultDict

import itertools
from collections import defaultdict
from datetime import datetime, timedelta
import click
import requests
from requests.exceptions import ConnectTimeout, MissingSchema
from sqlalchemy import func, and_

from scripts.db.database import session_handler
from scripts.db.models import Sample

STRAIN_LEVEL_AND_GLOBAL_CONTEXT_REPORT_HEADERS = ("Lineage", "Description", "Count", "Percent", "Global Count")
STRAIN_FIRST_SEEN_HEADERS = (
    "Pangolin Lineage",
    "Date First Seen",
    "MRN",
    "Sample ID",
    "Area",
    "Block",
    "Travel History",
    "Number of Subsequent occurrences",
)
STRAIN_PREVALENCE_HEADERS = (
    "Location",
    "Dominant Lineage",
    "Frequency",
    "Frequency Last Week",
    "Difference",
    "Sample Number",
)
STRAIN_PREVALENCE_WEEKS = 2
FLOAT_ROUNDING = 2
REPORT_RETURN = List[Tuple[Any, ...]]
GROUPING_DICT = DefaultDict[str, List[List[Any]]]
LINEAGES_COUNT_DICT = Dict[str, Dict[str, int]]


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
        content = (x.replace(f"{delimiter}{delimiter}", delimiter) for x in content)
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
) -> REPORT_RETURN:
    """
    Function to fetch report data from local database and external data sources.
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
    result: REPORT_RETURN = [STRAIN_LEVEL_AND_GLOBAL_CONTEXT_REPORT_HEADERS]
    notes = dict(notes)
    total = sum(x[1] for x in lineage_count)
    global_total = 0
    for lineage in lineage_count:
        name = lineage[0] or ""
        global_count = sum(1 for x in meta_data if x[-1] == name)
        global_total += global_count
        percentage = "{:.0%}".format(lineage[1] / total)
        result.append(
            (
                name,
                notes.get(name, "None available"),
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


def get_strain_first_seen_report_data() -> REPORT_RETURN:
    """
    Function to fetch report data from local database.
    :return List[Tuple[Any, ...]]: lines for the csv report containing columns corresponding
    to STRAIN_FIRST_SEEN_HEADERS
    """

    def cleanup_result(raw_data, counts_dict):
        res = []
        for line in raw_data:
            date_collected = line[1] or ""
            if isinstance(date_collected, datetime):
                date_collected = date_collected.date()
            res.append(
                (
                    line[0] or "",
                    date_collected,
                    line[2] or "",
                    line[3] or "",
                    line[4] or "",
                    line[5] or "",
                    line[6] or "",
                    counts_dict[line[0]] - 1,  # reducing by 1 for column "Number of Subsequent occurrences",
                )
            )
        return res

    result: REPORT_RETURN = [STRAIN_FIRST_SEEN_HEADERS]

    with session_handler() as session:
        counts = dict(session.query(Sample.pangolin_lineage, func.count("*")).group_by(Sample.pangolin_lineage).all())
        samples = (
            session.query(
                Sample.pangolin_lineage,
                Sample.date_collected,
                Sample.mrn,
                Sample.lab_id,
                Sample.area_name,
                Sample.block_number,
                Sample.travel_exposure,
            )
            .distinct(Sample.pangolin_lineage)
            .order_by(Sample.pangolin_lineage, Sample.date_collected)
            .all()
        )
        result.extend(cleanup_result(samples, counts))
    return result


def get_strain_prevalence_report_data() -> REPORT_RETURN:
    """
    Function to fetch report data from local database.
    :return List[Tuple[Any, ...]]: lines for the csv report containing columns corresponding
    to STRAIN_PREVALENCE_HEADERS
    """
    week_ago: datetime = datetime.now() - timedelta(weeks=1)
    result: REPORT_RETURN = [STRAIN_PREVALENCE_HEADERS]

    def get_frequency(count, total):
        if total:
            return round(count / total, FLOAT_ROUNDING)
        return 0

    def process_lineages(lineage_records: List[List[Any]]) -> LINEAGES_COUNT_DICT:
        """
        group records by lineages and count how many are
        from the last week, older than one week, and total
        """
        totals = {}
        groups: GROUPING_DICT = defaultdict(list)

        for lineage, lines in itertools.groupby(lineage_records, lambda x: x[0]):
            groups[lineage].extend(list(line) for line in lines)

        for lineage, group in groups.items():
            totals[lineage] = {
                "this_week": sum(1 for x in group if x[3] > week_ago),
                "last_week": sum(1 for x in group if x[3] <= week_ago),
                "count": len(group),
            }
        return totals

    def get_prevalent_lineage(lineages, total_this_week, total_last_week):
        """
        determines prevalent strain commencing STRAIN_PREVALENCE_WEEKS prior to today
        """
        prevalent = ""
        prevalent_count = {
            "last_week": 0,
            "count": 0,
        }
        last_week_frequency = 0
        for lineage, counts in lineages.items():
            if counts["count"] < prevalent_count["count"]:
                continue

            last_week_frequency = get_frequency(counts["last_week"], total_last_week)
            if counts["count"] == prevalent_count["count"] and last_week_frequency > get_frequency(
                prevalent_count["last_week"], total_last_week
            ):
                continue

            prevalent_count = counts
            prevalent = lineage

        count_this_week = prevalent_count.get("this_week", 0)
        frequency = get_frequency(count_this_week, total_this_week)
        return prevalent, count_this_week, frequency, last_week_frequency

    def count_by_column(records: List[Any], column: int = None) -> None:
        """
        groups `records` data by specified `column`, determines `prevelant_lineage` for each group,
        and calculates frequencies. If `column` is `None`, calculates entire country from all data given
        """
        groups: GROUPING_DICT = defaultdict(list)

        if column is None:
            groups["Bahrain"] = records
        else:
            for location, values in itertools.groupby(records, key=lambda x: x[column]):
                groups[location].extend(list(values))

        for location, group in groups.items():
            total_this_week: int = sum(1 for x in group if x[3] > week_ago)
            total_last_week: int = sum(1 for x in group if x[3] <= week_ago)
            lineages: LINEAGES_COUNT_DICT = process_lineages(group)

            prevalent, count, frequency, last_week_frequency = get_prevalent_lineage(
                lineages, total_this_week, total_last_week
            )
            if not prevalent:
                continue

            result.append(
                (
                    location,
                    prevalent,
                    frequency,
                    last_week_frequency,
                    round(frequency - last_week_frequency, FLOAT_ROUNDING),
                    count,
                )
            )

    with session_handler() as session:
        samples = (
            session.query(Sample.pangolin_lineage, Sample.area_name, Sample.governorate_name, Sample.date_collected)
            .filter(
                and_(
                    Sample.date_collected >= func.now() - timedelta(weeks=STRAIN_PREVALENCE_WEEKS),
                    Sample.pangolin_lineage.isnot(None),
                )
            )
            .all()
        )

    # area
    count_by_column(samples, 1)
    # governorate
    count_by_column(samples, 2)
    # entire country
    count_by_column(samples)

    return result


reports: Dict[str, Callable[..., REPORT_RETURN]] = {
    "strain_level_and_global_context": get_strain_level_and_global_context_report_data,
    "strain_first_seen": get_strain_first_seen_report_data,
    "strain_prevalence": get_strain_prevalence_report_data,
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
    type=click.Choice(["strain_level_and_global_context", "strain_first_seen", "strain_prevalence"]),
    required=True,
    help="name of the report to be outputted",
)
@click.option(
    "--lineage-notes-url",
    type=click.STRING,
    help="url to latest pangoLEARN lineage notes",
)
@click.option(
    "--metadata-url",
    type=click.STRING,
    help="url to latest pangoLEARN metadata",
)
@click.option(
    "--pangolearn-dir",
    type=click.Path(dir_okay=True, readable=True),
    help="fallback if latest data fetch fails, will look for latest local data by filename from url",
)
def generate_report(output: str, report: str, **kwargs) -> None:
    """
    Generate a csv report from the database and external sources
    """
    kwargs = {key: value for key, value in kwargs.items() if value}
    samples = reports[report](**kwargs)

    with open(output, "w") as out_file:
        csv_writer = csv.writer(out_file)
        csv_writer.writerows(samples)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_report()
