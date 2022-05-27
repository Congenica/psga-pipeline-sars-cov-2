import pytest
import structlog

from scripts.util.logging import LOG_LEVELS, get_structlog_logger
from tests.util.test_notification import load_log_file_to_dict


@pytest.mark.parametrize(
    "message,sample",
    [
        (
            "Here is a log message",
            "silly_sample_for_reusing_the_function",
        ),
    ],
)
@pytest.mark.parametrize(
    "level",
    [level.lower() for level in LOG_LEVELS],
)
def test_logger(tmp_path, message, level, sample):
    log_file = tmp_path / "messages.log"
    assert not log_file.is_file()

    structlog.reset_defaults()
    logger_level = getattr(get_structlog_logger(log_file=f"{log_file}", log_level="DEBUG"), level)
    logger_level(message, sample=sample)

    # test the file handler
    assert log_file.is_file()
    log_dict = load_log_file_to_dict(log_file, "sample")
    log_sample = log_dict[sample]
    assert sample == log_sample["sample"]
    assert level == log_sample["level"]
    assert message == log_sample["message"]


def test_logger_splunk_handler(tmp_path, setup_splunk):
    message = "a message to be logged to Splunk"

    log_file = tmp_path / "messages.log"
    assert not log_file.is_file()

    structlog.reset_defaults()
    logger = get_structlog_logger(log_file=f"{log_file}")

    logger.info(message)
