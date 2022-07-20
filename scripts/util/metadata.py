import csv
from uuid import UUID
from typing import List
from pathlib import Path
from dataclasses import dataclass, field
import click

from scripts.util.notifications import Event, Notification
from scripts.validation.check_csv_columns import check_csv_columns

EXPECTED_HEADERS = {"sample_id", "file_1", "file_2"}


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

    if sample_with_two_reads:
        if not row["file_2"]:
            errs.append(f"file_2 for {sample_id} not available")

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


def inspect_metadata_file(
    metadata_path: Path,
    samples_with_two_reads: bool,
) -> ProcessedSamples:
    """
    Inspect the metadata file
    """
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
                processed_samples.valid.append(sample_name)

            except ValueError as e:
                sample_name = row["sample_id"]
                click.echo(f"Invalid row for sample {sample_name}:\n{e}", err=True)
                processed_samples.invalid.append(sample_name)

    return processed_samples
