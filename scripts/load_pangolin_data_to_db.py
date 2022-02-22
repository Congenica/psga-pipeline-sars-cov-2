from typing import Dict
import csv

import click
from click import ClickException
from sqlalchemy.orm import scoped_session

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample, PangolinStatus

EXPECTED_HEADERS = {
    "taxon",
    "status",
    "lineage",
    "conflict",
    "ambiguity_score",
    "scorpio_call",
    "scorpio_support",
    "scorpio_conflict",
    "note",
    "pangolin_version",
    "pangoLEARN_version",
    "pango_version",
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

    pangolin_status = PangolinStatus[sample_from_csv["status"]]

    if pangolin_status == PangolinStatus.passed_qc:
        sample.pangolin_lineage = sample_from_csv["lineage"]

    sample.pangolin_conflict = sample_from_csv["conflict"] if sample_from_csv["conflict"] else None
    sample.pangolin_ambiguity_score = sample_from_csv["ambiguity_score"] if sample_from_csv["ambiguity_score"] else None
    sample.pangolin_status = pangolin_status

    sample.scorpio_call = sample_from_csv["scorpio_call"] if sample_from_csv["scorpio_call"] else None
    sample.scorpio_support = sample_from_csv["scorpio_support"] if sample_from_csv["scorpio_support"] else None
    sample.scorpio_conflict = sample_from_csv["scorpio_conflict"] if sample_from_csv["scorpio_conflict"] else None
    sample.note = sample_from_csv["note"] if sample_from_csv["note"] else None


@click.command()
@click.option(
    "--pangolin-lineage-report-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="A CSV report combining all samples' Pangolin pipeline output lineage reports",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
def load_pangolin_data(pangolin_lineage_report_file: str, analysis_run_name: str) -> None:
    """
    Load Pangolin lineage report for a certain sample to the database
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

        with session_handler() as session:

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
                analysis_run.pangolin_version = record["pangolin_version"] if record["pangolin_version"] else None
                analysis_run.pangolearn_version = record["pangoLEARN_version"] if record["pangoLEARN_version"] else None
                analysis_run.pango_version = record["pango_version"] if record["pango_version"] else None


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_pangolin_data()
