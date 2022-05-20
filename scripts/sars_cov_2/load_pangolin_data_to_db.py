from typing import Dict, List
from pathlib import Path
import csv
from distutils.util import strtobool

import click
from click import ClickException
from sqlalchemy.orm import scoped_session

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample, PangolinStatus
from scripts.util.notifications import Notification

EXPECTED_HEADERS = {
    "taxon",
    "lineage",
    "conflict",
    "ambiguity_score",
    "scorpio_call",
    "scorpio_support",
    "scorpio_conflict",
    "scorpio_notes",
    "version",
    "pangolin_version",
    "scorpio_version",
    "constellation_version",
    "is_designated",
    "qc_status",
    "qc_notes",
    "note",
}


def load_pangolin_sample(session: scoped_session, analysis_run_name: str, sample_from_csv: Dict):

    sample_name = sample_from_csv["taxon"]

    sample = (
        session.query(Sample)
        .join(AnalysisRun)
        .filter(
            Sample.sample_name == sample_name,
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )

    if not sample:
        # This should never happen as samples are expected to be loaded as part of the metadata loading
        # if this happens, there is an error in the pipeline
        raise ClickException(f"Sample name: {sample_name} for analysis run: {analysis_run_name} was not found")

    pangolin_status = PangolinStatus[sample_from_csv["qc_status"].upper()]

    if pangolin_status == PangolinStatus.PASS:
        sample.pangolin_lineage = sample_from_csv["lineage"]

    sample.pangolin_conflict = sample_from_csv["conflict"] if sample_from_csv["conflict"] else None
    sample.pangolin_ambiguity_score = sample_from_csv["ambiguity_score"] if sample_from_csv["ambiguity_score"] else None
    sample.pangolin_status = pangolin_status

    sample.scorpio_call = sample_from_csv["scorpio_call"] if sample_from_csv["scorpio_call"] else None
    sample.scorpio_support = sample_from_csv["scorpio_support"] if sample_from_csv["scorpio_support"] else None
    sample.scorpio_conflict = sample_from_csv["scorpio_conflict"] if sample_from_csv["scorpio_conflict"] else None
    sample.scorpio_notes = sample_from_csv["scorpio_notes"] if sample_from_csv["scorpio_notes"] else None
    sample.is_designated = (
        bool(strtobool(sample_from_csv["is_designated"])) if sample_from_csv["is_designated"] else None
    )
    sample.qc_notes = sample_from_csv["qc_notes"] if sample_from_csv["qc_notes"] else None
    sample.note = sample_from_csv["note"] if sample_from_csv["note"] else None


def load_pangolin_lineages(session: scoped_session, pangolin_lineage_report_file: Path, analysis_run_name: str) -> None:
    """
    Load pangolin sample data and pangolin version to the database
    """

    with open(pangolin_lineage_report_file) as csv_file:
        sample_from_csv_reader = csv.DictReader(csv_file)

        headers = set(sample_from_csv_reader.fieldnames)
        if not EXPECTED_HEADERS.issubset(headers):
            err = (
                "Unexpected CSV headers, got:\n"
                + ", ".join(headers)
                + "\n, but expect at least \n"
                + ", ".join(EXPECTED_HEADERS)
            )
            raise ClickException(err)

        # update pangolin-related data in analysis_run table
        analysis_run = (
            session.query(AnalysisRun).filter(AnalysisRun.analysis_run_name == analysis_run_name).one_or_none()
        )
        if not analysis_run:
            raise ClickException(f"Analysis run {analysis_run_name} was not found in the database")

        # grab the first record in order to get the pango- versions
        record: Dict = None

        for idx, sample_from_csv in enumerate(sample_from_csv_reader):
            if idx == 0:
                record = sample_from_csv
            load_pangolin_sample(session, analysis_run_name, sample_from_csv)

        # these are the same for all samples
        if record:
            analysis_run.constellation_version = (
                record["constellation_version"] if record["constellation_version"] else None
            )
            analysis_run.pangolin_version = record["pangolin_version"] if record["pangolin_version"] else None
            analysis_run.pangolin_data_version = record["version"] if record["version"] else None
            analysis_run.scorpio_version = record["scorpio_version"] if record["scorpio_version"] else None


def get_analysis_run_samples(session: scoped_session, analysis_run_name: str) -> List[Sample]:
    """
    Return all the samples for this analysis run
    """
    samples = session.query(Sample).join(AnalysisRun).filter(AnalysisRun.analysis_run_name == analysis_run_name).all()
    return samples


def generate_notifications(
    session: scoped_session,
    analysis_run_name: str,
    no_qc_path: Path,
    failed_qc_path: Path,
    passed_qc_path: Path,
) -> None:
    """
    Generate and publish the notifications for pangolin
    """
    samples = get_analysis_run_samples(session, analysis_run_name)

    notifications = Notification(
        events={
            "no_qc": {
                "path": no_qc_path,
                "level": "ERROR",
                "message": "pangolin did not terminate successfully for this sample",
                "samples": [s.sample_name for s in samples if s.pangolin_status == PangolinStatus.UNKNOWN],
            },
            "failed_qc": {
                "path": failed_qc_path,
                "level": "WARNING",
                "event": "pangolin failed",
                "samples": [s.sample_name for s in samples if s.pangolin_status == PangolinStatus.FAIL],
            },
            "passed_qc": {
                "path": passed_qc_path,
                "level": "INFO",
                "event": "pangolin qc passed",
                "samples": [s.sample_name for s in samples if s.pangolin_status == PangolinStatus.PASS],
            },
        }
    )

    notifications.publish()


@click.command()
@click.option(
    "--pangolin-lineage-report-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="A CSV report combining all samples' Pangolin pipeline output lineage reports",
)
@click.option(
    "--samples-with-unknown-pangolin-status-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the names of samples for this analysis run which have UNKNOWN pangolin status",
)
@click.option(
    "--samples-with-failed-pangolin-status-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the names of samples for this analysis run which have FAIL pangolin status",
)
@click.option(
    "--samples-with-passed-pangolin-status-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the names of samples for this analysis run which have PASS pangolin status",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
def load_pangolin_data(
    pangolin_lineage_report_file: str,
    samples_with_unknown_pangolin_status_file: str,
    samples_with_failed_pangolin_status_file: str,
    samples_with_passed_pangolin_status_file: str,
    analysis_run_name: str,
) -> None:
    """
    Load Pangolin lineage report for a certain sample to the database
    """
    with session_handler() as session:
        load_pangolin_lineages(session, Path(pangolin_lineage_report_file), analysis_run_name)

        generate_notifications(
            session,
            analysis_run_name,
            Path(samples_with_unknown_pangolin_status_file),
            Path(samples_with_failed_pangolin_status_file),
            Path(samples_with_passed_pangolin_status_file),
        )


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_pangolin_data()
