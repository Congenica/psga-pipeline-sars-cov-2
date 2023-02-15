from dataclasses import dataclass, field

from scripts.util.notifications import Event, Notification

SAMPLE_ID = "SAMPLE_ID"
SEQ_FILE_1 = "SEQ_FILE_1"
SEQ_FILE_2 = "SEQ_FILE_2"

EXPECTED_HEADERS = {SAMPLE_ID, SEQ_FILE_1, SEQ_FILE_2}
ILLUMINA = "illumina"
ONT = "ont"
UNKNOWN = "unknown"


@dataclass
class ProcessedSamples:
    valid: list[str] = field(metadata={"required": True}, default_factory=list)
    invalid: list[str] = field(metadata={"required": True}, default_factory=list)


def normalise_row(row):
    """
    Strip leading and trailing spaces from everything
    """
    for f in EXPECTED_HEADERS:
        row[f] = row[f].lstrip().rstrip() if row[f] is not None else ""
    return row


def generate_notifications(
    analysis_run_name: str,
    valid_samples: list[str],
    invalid_samples: list[str],
) -> None:
    """
    Generate and publish the notifications for ncov
    """
    notifications = Notification(
        events={
            "failed_qc": Event(
                analysis_run=analysis_run_name,
                level="ERROR",
                message="metadata validation failed",
                samples=invalid_samples,
            ),
            "passed_qc": Event(
                analysis_run=analysis_run_name,
                level="INFO",
                message="metadata validation passed",
                samples=valid_samples,
            ),
        }
    )

    notifications.publish()
