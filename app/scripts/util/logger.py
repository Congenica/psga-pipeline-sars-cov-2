import json

# pylint: disable=import-self
import logging
import structlog

CRITICAL = "CRITICAL"
ERROR = "ERROR"
WARNING = "WARNING"
INFO = "INFO"
DEBUG = "DEBUG"

LOG_LEVELS = {CRITICAL, ERROR, WARNING, INFO, DEBUG}


def key_ordering_serialiser(data: dict, **kwargs) -> str:
    """
    Reorder the dictionary to make the log lines a bit easier to scan by eye.
    :param data: dictionary of data to log
    :return: data converted to json string after reordering.
    """
    ordered_keys = ["timestamp", "level", "logger"]
    ordered_data: dict = {}
    # Add the keys in the order specified by the list.
    for key in ordered_keys:
        ordered_data[key] = data.pop(key)

    # Add the remaining keys in their existing order.
    ordered_data.update(data)

    # finally rename "event" to "message" as that is more common
    if "event" in ordered_data:
        ordered_data["message"] = ordered_data.pop("event")

    return json.dumps(ordered_data, **kwargs)


def create_structlog_processors() -> list:
    """
    Create the list of log even processors for structlog.
    The processors are applied in order to a logged event, and add useful info such as
    timestamps and logging levels. The last processor is the structlog renderer, which will turn
    all logged data into a JSON object.
    :return: A list of structlog processors.
    """
    processors = [
        structlog.stdlib.filter_by_level,
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.stdlib.add_log_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.PositionalArgumentsFormatter(),
        # Allows easy passing and formatting of stack traces in log messages.
        structlog.processors.StackInfoRenderer(),
        # Replace exc_info with nicely formatted exception.
        structlog.processors.format_exc_info,
        # JSON rendering
        structlog.processors.JSONRenderer(serializer=key_ordering_serialiser),
    ]
    return processors


def get_structlog_logger(
    logger_name: str = "psga", log_file: str = "psga.log", log_level: str = "INFO"
) -> structlog.wrap_logger:
    """
    This method sets up structlog. It must be called before anything else attempts to
    initialise or use the stdlib logger.
    :param logger_name: The name of the logger.
    :param log_file: The name of the file to log.
    :param log_level: The log level for the logging.
    :return: Logging configuration dictionary.
    """
    if structlog.is_configured():
        # avoid: RuntimeWarning `repeated configuration attempted`
        return structlog.get_logger(logger_name)

    # Setup the logger with its handlers first, as structlog wraps it
    logger = logging.getLogger(logger_name)  # type: ignore

    logger.setLevel(log_level)
    formatter = logging.Formatter("%(message)s")  # type: ignore

    # stream handler
    stream_handler = logging.StreamHandler()  # type: ignore
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    # file handler
    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")  # type: ignore
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Setup structlog
    processors = create_structlog_processors()
    structlog.configure_once(
        processors=processors,
        # Use thread local storage so the logging context is global to the request thread.
        context_class=structlog.threadlocal.wrap_dict(dict),
        # Interface with the standard library logger which we configured above.
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )

    struct_logger = structlog.wrap_logger(logger)

    return struct_logger
