from typing import Dict
import csv

import click
from sqlalchemy.orm import scoped_session, joinedload

from db.database import session_handler
from db.models import Sample, SampleQC, PangolinStatus


def load_data_from_csv(session: scoped_session, sample_name: str, sample_from_csv: Dict):
    sample = (
        session.query(Sample).filter_by(sample_name=sample_name).options(joinedload(Sample.sample_qc)).one_or_none()
    )

    if not sample:
        sample = Sample(sample_name=sample_name)
        session.add(sample)
    if not sample.sample_qc:
        sample.sample_qc = SampleQC()

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

    sample.sample_qc.pangolin_version = sample_from_csv["pangolin_version"]
    sample.sample_qc.pangolearn_version = sample_from_csv["pangoLEARN_version"]
    sample.sample_qc.pango_version = sample_from_csv["pango_version"]


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
def load_pangolin_data(pangolin_lineage_report_file: str, sample_name: str) -> None:
    """
    Load Pangolin lineage report for a certain sample to the database
    """
    with open(pangolin_lineage_report_file) as csv_file:
        sample_from_csv_reader = csv.DictReader(csv_file)
        with session_handler() as session:
            # this file only has one data row
            for sample_from_csv in sample_from_csv_reader:
                load_data_from_csv(session, sample_name, sample_from_csv)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_pangolin_data()
