import os
import csv
from pathlib import Path
from typing import Dict, List
import click
from click import ClickException
from sqlalchemy.orm import scoped_session, joinedload

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample, SampleQC
from scripts.util.notifications import Notification
from scripts.validation.check_csv_columns import check_csv_columns

EXPECTED_HEADERS = {
    "sample_name",
    "pct_N_bases",
    "pct_covered_bases",
    "longest_no_N_run",
    "num_aligned_reads",
    "qc_pass",
}


def load_ncov_sample(
    session: scoped_session, analysis_run_name: str, sample_from_csv: Dict, ncov_qc_depth_directory: str
):
    """
    Load the ncov results of a sample to the database
    """
    sample_name = sample_from_csv["sample_name"]
    sample_qc_depth_file_path = os.path.join(ncov_qc_depth_directory, f"{sample_name}.depth.png")
    if not os.path.isfile(sample_qc_depth_file_path):
        raise ValueError(f"File {sample_qc_depth_file_path}, required to submit sample {sample_name}, does not exist")

    sample = (
        session.query(Sample)
        .join(AnalysisRun)
        .filter(
            Sample.sample_name == sample_name,
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .options(joinedload(Sample.sample_qc))
        .one_or_none()
    )

    if not sample:
        raise ClickException(f"Sample name: {sample_name} was not found")

    if not sample.sample_qc:
        sample.sample_qc = SampleQC()

    sample_qc = sample.sample_qc
    sample_qc.pct_n_bases = sample_from_csv["pct_N_bases"]
    sample_qc.pct_covered_bases = sample_from_csv["pct_covered_bases"]
    sample_qc.longest_no_n_run = sample_from_csv["longest_no_N_run"]
    sample_qc.num_aligned_reads = sample_from_csv["num_aligned_reads"]
    sample_qc.qc_pass = sample_from_csv["qc_pass"].lower() == "true"
    with open(sample_qc_depth_file_path, "rb") as f:
        sample_qc.qc_plot = bytearray(f.read())


def get_analysis_run_samples(session: scoped_session, analysis_run_name: str) -> List[Sample]:
    """
    Return all the samples for this analysis run
    """
    samples = (
        session.query(Sample)
        .join(AnalysisRun)
        .filter(AnalysisRun.analysis_run_name == analysis_run_name)
        .options(joinedload(Sample.sample_qc))
        .all()
    )
    return samples


def generate_notifications(
    session: scoped_session,
    analysis_run_name: str,
    no_qc_path: Path,
    failed_qc_path: Path,
    passed_qc_path: Path,
) -> None:
    """
    Generate and publish the notifications for ncov
    """
    samples = get_analysis_run_samples(session, analysis_run_name)

    notifications = Notification(
        events={
            "no_qc": {
                "path": no_qc_path,
                "level": "ERROR",
                "message": "ncov did not terminate successfully for this sample",
                "samples": [s.sample_name for s in samples if not s.sample_qc],
            },
            "failed_qc": {
                "path": failed_qc_path,
                "level": "WARNING",
                "event": "ncov qc failed",
                "samples": [s.sample_name for s in samples if s.sample_qc and not s.sample_qc.qc_pass],
            },
            "passed_qc": {
                "path": passed_qc_path,
                "level": "INFO",
                "event": "ncov qc passed",
                "samples": [s.sample_name for s in samples if s.sample_qc and s.sample_qc.qc_pass],
            },
        }
    )

    notifications.publish()


@click.command()
@click.option(
    "--ncov-qc-csv-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="ncov pipeline resulting qc csv file",
)
@click.option(
    "--ncov-qc-depth-directory",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    required=True,
    help="directory, containing qc depth files, following the pattern {sample_name}.depth.png",
)
@click.option(
    "--samples-without-ncov-qc-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the names of samples for this analysis run which do not have any sample_qc record",
)
@click.option(
    "--samples-with-failed-ncov-qc-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the names of samples for this analysis run which have failed ncov qc",
)
@click.option(
    "--samples-with-passed-ncov-qc-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file storing the names of samples for this analysis run which have passed ncov qc",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
def load_ncov_data(
    ncov_qc_csv_file: str,
    ncov_qc_depth_directory: str,
    samples_without_ncov_qc_file: str,
    samples_with_failed_ncov_qc_file: str,
    samples_with_passed_ncov_qc_file: str,
    analysis_run_name: str,
) -> None:
    """
    Submit samples QC to the database, generated by ncov pipeline
    """
    with session_handler() as session:
        with open(ncov_qc_csv_file) as csv_file:
            sample_from_csv_reader = csv.DictReader(csv_file)
            check_csv_columns(
                set(sample_from_csv_reader.fieldnames),
                EXPECTED_HEADERS,
            )
            for sample_from_csv in sample_from_csv_reader:
                load_ncov_sample(session, analysis_run_name, sample_from_csv, ncov_qc_depth_directory)

        generate_notifications(
            session,
            analysis_run_name,
            Path(samples_without_ncov_qc_file),
            Path(samples_with_failed_ncov_qc_file),
            Path(samples_with_passed_ncov_qc_file),
        )


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_ncov_data()
