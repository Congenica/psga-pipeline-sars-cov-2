from uuid import UUID
from typing import List
from pathlib import Path

from scripts.util.notifications import Notification

METADATA_FILE_EXPECTED_HEADERS = {"sample_id", "file_1", "file_2", "md5_1", "md5_2"}


def is_valid_uuid(uuid_str: str) -> bool:
    try:
        UUID(uuid_str)
        return True
    except ValueError:
        return False


def generate_notifications(
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
            "failed_qc": {
                "path": failed_qc_path,
                "level": "ERROR",
                "event": "metadata validation failed",
                "samples": invalid_samples,
            },
            "passed_qc": {
                "path": passed_qc_path,
                "level": "INFO",
                "event": "metadata validation passed",
                "samples": valid_samples,
            },
        }
    )

    notifications.publish()
