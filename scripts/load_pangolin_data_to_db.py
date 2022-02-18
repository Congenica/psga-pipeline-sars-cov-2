from typing import Dict
import csv

import click
from click import ClickException
from sqlalchemy.orm import scoped_session

from db.database import session_handler
from db.models import AnalysisRun, Sample, PangolinStatus


def load_data_from_csv(session: scoped_session, analysis_run_name: str, sample_name: str, sample_from_csv: Dict):
    sample = (
        session.query(Sample)
        .filter(
            Sample.sample_name == sample_name,
        )
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )

    sample_analysis_run = (
        session.query(AnalysisRun).filter(AnalysisRun.analysis_run_id == sample.analysis_run_id).one_or_none()
    )

    if not sample:
        raise ClickException(f"Sample name: {sample_name} was not found")

    if not sample_analysis_run:
        raise ClickException(f"No analysis run was found for sample name: {sample_name}")

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

    # This needs refactory. We should run all pangolin data in one single process
    # and this part should be done one time only
    sample_analysis_run.pangolin_version = (
        sample_from_csv["pangolin_version"] if sample_from_csv["pangolin_version"] else None
    )
    sample_analysis_run.pangolearn_version = (
        sample_from_csv["pangoLEARN_version"] if sample_from_csv["pangoLEARN_version"] else None
    )
    sample_analysis_run.pango_version = sample_from_csv["pango_version"] if sample_from_csv["pango_version"] else None


@click.command()
@click.option(
    "--pangolin-lineage-report-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Pangolin pipeline output lineage report with pattern {sample_name}_lineage_report.csv",
)
@click.option(
    "--sample-name",
    type=str,
    required=True,
    help="Lab sample identifier",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
def load_pangolin_data(pangolin_lineage_report_file: str, sample_name: str, analysis_run_name: str) -> None:
    """
    Load Pangolin lineage report for a certain sample to the database
    """
    with open(pangolin_lineage_report_file) as csv_file:
        sample_from_csv_reader = csv.DictReader(csv_file)
        with session_handler() as session:
            # this file only has one data row
            for sample_from_csv in sample_from_csv_reader:
                load_data_from_csv(session, analysis_run_name, sample_name, sample_from_csv)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_pangolin_data()
