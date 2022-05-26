import csv
from uuid import UUID
from typing import List
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field
from sqlalchemy.orm import scoped_session
import click

from scripts.db.models import AnalysisRun, Sample
from scripts.db.queries import get_analysis_run_sample
from scripts.util.notifications import Event, Notification
from scripts.validation.check_csv_columns import check_csv_columns

EXPECTED_HEADERS = {"sample_id", "file_1", "file_2", "md5_1", "md5_2"}


@dataclass
class ProcessedSamples:
    valid: List[str] = field(metadata={"required": True}, default_factory=list)
    invalid: List[str] = field(metadata={"required": True}, default_factory=list)


def is_valid_uuid(uuid_str: str) -> bool:
    try:
        UUID(uuid_str)
        return True
    except ValueError:
        return False


def validate_and_normalise_row(row, sample_with_two_reads: bool):

    # strip leading and trailing spaces from everything
    for f in EXPECTED_HEADERS:
        row[f] = row[f].lstrip().rstrip() if row[f] is not None else ""

    errs = []

    # add validations here
    sample_id = row["sample_id"]
    if not sample_id:
        errs.append("sample_id not available")
    elif not is_valid_uuid(sample_id):
        errs.append(f'sample_id "{sample_id}" is not a UUID')

    if not row["file_1"]:
        errs.append(f"file_1 for {sample_id} not available")
    if not row["md5_1"]:
        errs.append(f"md5_1 for {sample_id} not available")

    if sample_with_two_reads:
        if not row["file_2"]:
            errs.append(f"file_2 for {sample_id} not available")
        if not row["md5_2"]:
            errs.append(f"md5_2 for {sample_id} not available")

    if errs:
        raise ValueError("\n".join(errs))

    return row


def generate_notifications(
    analysis_run_name: str,
    valid_samples: List[str],
    passed_qc_path: Path,
    invalid_samples: List[str],
    failed_qc_path: Path,
) -> None:
    """
    Generate and publish the notifications for ncov
    """
    notifications = Notification(
        events={
            "failed_qc": Event(
                analysis_run=analysis_run_name,
                path=failed_qc_path,
                level="ERROR",
                message="metadata validation failed",
                samples=invalid_samples,
            ),
            "passed_qc": Event(
                analysis_run=analysis_run_name,
                path=passed_qc_path,
                level="INFO",
                message="metadata validation passed",
                samples=valid_samples,
            ),
        }
    )

    notifications.publish()


def link_sample_to_analysis_run(
    session: scoped_session,
    analysis_run: AnalysisRun,
    sample_name: str,
    load_missing_samples: bool = False,
) -> None:
    """
    Link a sample to the analysis run.
    Raise a ValueError if the sample does not exist in the database,
    unless load_missing_samples is True.
    """

    sample = get_analysis_run_sample(session, analysis_run.analysis_run_name, sample_name)

    if sample:
        sample.analysis_run_id = analysis_run.analysis_run_id

    elif load_missing_samples:
        # the pipeline metadata does not contain the date_collected field
        yesterday = datetime.strftime(datetime.now(), "%Y-%m-%d")
        sample = Sample(
            analysis_run_id=analysis_run.analysis_run_id,
            sample_name=sample_name,
            date_collected=yesterday,
            metadata_loaded=True,
        )
        session.add(sample)

    else:
        raise ValueError(f"Sample {sample_name} not found in the database, but listed in pipeline metadata")

    session.commit()


def inspect_metadata_file(
    session: scoped_session,
    metadata_path: Path,
    analysis_run: AnalysisRun,
    samples_with_two_reads: bool,
    load_missing_samples: bool,
) -> ProcessedSamples:
    """
    Load the metadata and return a csv.DictReader object for this metadata.
    Raise a ValueError exception if the metadata contains less columns than
    what expected.
    """

    if not analysis_run:
        raise ValueError("Cannot retrieve samples for a non existing analysis run.")

    processed_samples = ProcessedSamples(
        valid=[],
        invalid=[],
    )

    with open(metadata_path, "r") as metadata_fd:
        reader = csv.DictReader(metadata_fd, delimiter=",")

        check_csv_columns(set(reader.fieldnames), EXPECTED_HEADERS)

        for row in reader:
            try:
                row = validate_and_normalise_row(row, samples_with_two_reads)
                sample_name = row["sample_id"]
                link_sample_to_analysis_run(session, analysis_run, sample_name, load_missing_samples)

                processed_samples.valid.append(sample_name)

            except ValueError as e:
                sample_name = row["sample_id"]
                click.echo(f"Invalid row for sample {sample_name}:\n{e}", err=True)
                processed_samples.invalid.append(sample_name)

    return processed_samples
