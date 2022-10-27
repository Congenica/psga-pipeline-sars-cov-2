from pathlib import Path
from typing import Dict
import json
import pytest
import structlog

from scripts.util.notifications import Event, Notification
from scripts.util.logging import get_structlog_logger, ERROR, INFO


def load_log_file_to_dict(log_file: Path, key: str) -> Dict:
    log_dict = {}
    with open(log_file, "r") as fd:
        for line in fd.read().splitlines():
            line_dict = json.loads(line)
            log_dict[line_dict[key]] = line_dict
    return log_dict


@pytest.mark.parametrize(
    "events",
    [
        {},
        {
            "evt1": Event(
                analysis_run="analysis_run1",
                level=ERROR,
                message="error message",
                samples=["a", "b"],
            ),
        },
        {
            "evt1": Event(
                analysis_run="analysis_run1",
                level=INFO,
                message="info message",
                samples=["x", "y", "z"],
            ),
        },
        {
            "evt1": Event(
                analysis_run="analysis_run1",
                level=ERROR,
                message="error message",
                samples=["a", "b"],
            ),
            "evt2": Event(
                analysis_run="analysis_run1",
                level=INFO,
                message="info message",
                samples=["x", "y", "z"],
            ),
        },
        {
            "evt1": Event(
                analysis_run="analysis_run1",
                level=ERROR,
                message="error message",
                samples=["a", "b"],
            ),
            "evt2": Event(
                analysis_run="analysis_run2",
                level=INFO,
                message="info message",
                samples=["x", "y", "z"],
            ),
        },
    ],
)
def test_notification(tmp_path, events):

    log_file = tmp_path / "messages.log"
    assert not log_file.is_file()

    structlog.reset_defaults()
    get_structlog_logger(log_file=f"{log_file}")

    notifications = Notification(events=events)
    notifications.publish()

    assert log_file.is_file()
    log_dict = load_log_file_to_dict(log_file, "sample")

    for evt in events.values():
        # assert log file
        for sample in evt.samples:
            log_sample = log_dict[sample]
            assert evt.level == log_sample["level"].upper()
            assert evt.message == log_sample["message"]
            assert evt.analysis_run == log_sample["analysis_run"]
            assert sample == log_sample["sample"]
