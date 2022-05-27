import os
import json
from typing import Dict, List

# pylint: disable=import-self
import logging
import structlog
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


LOG_LEVELS = {"CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"}


# pylint: disable=no-member
class HECHandler(logging.Handler):  # type: ignore
    """
    HEC logging handler.
    Posts logs to HEC endpoint (currently Splunk).
    Retries failed POST requests up to 3 times, but otherwise
    has no error handling: if log messages cannot be sent to
    the HEC an error is raised.

    :param endpoint: HEC endpoint
    :param token: authorisation token for endpoint
    :param hostname: identifier for the entity posting to the HEC (e.g. active pod name)
    :param source: identifier for source of logs (e.g. cluster)
    :param sourcetype: identifier for application sending logs (e.g. name of application)
    :param backoff_factor: the requests lib backoff factor
    :param total_retries: number of retry attempts on a failed/erroring connection
    """

    def __init__(
        self,
        endpoint: str,
        token: str,
        hostname: str,
        source: str,
        sourcetype: str,
        backoff_factor: float = 0.5,
        total_retries: int = 3,
    ) -> None:

        super().__init__()

        self.endpoint = endpoint
        self.token = token

        self.hostname = hostname
        self.source = source
        self.sourcetype = sourcetype

        self.session = requests.Session()
        self.session.headers.update(self.headers)

        retries = Retry(
            total=total_retries,
            backoff_factor=backoff_factor,
            status_forcelist=(500, 502, 503, 504),
            method_whitelist=["POST"],
        )

        self.session.mount(self.endpoint, HTTPAdapter(max_retries=retries))

    @property
    def headers(self) -> Dict[str, str]:
        """
        Authorisation (Splunk specific).
        """
        return {"Authorization": f"Splunk {self.token}"}

    def emit(self, record) -> None:
        """
        Send a log record to the HEC.
        """
        record = self.format(record)

        event = json.loads(record)
        payload = {"hostname": self.hostname, "source": self.source, "sourcetype": self.sourcetype, "event": event}

        r = self.session.post(self.endpoint, json=payload)
        r.raise_for_status()


def key_ordering_serialiser(data: Dict, **kwargs) -> str:
    """
    Reorder the dictionary to make the log lines a bit easier to scan by eye.
    :param data: Dictionary of data to log
    :return: data converted to json string after reordering.
    """
    ordered_keys = ["timestamp", "level", "logger"]
    ordered_data: Dict = {}
    # Add the keys in the order specified by the list.
    for key in ordered_keys:
        ordered_data[key] = data.pop(key)

    # Add the remaining keys in their existing order.
    ordered_data.update(data)

    # finally rename "event" to "message" as that is more common
    if "event" in ordered_data:
        ordered_data["message"] = ordered_data.pop("event")

    return json.dumps(ordered_data, **kwargs)


def create_structlog_processors() -> List:
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

    # hec handler
    if os.environ.get("HEC_ENABLED_CHILD") == "true":
        hostname = os.environ.get("HOSTNAME")
        if hostname is None:
            raise Exception("Cannot set logger hostname: HOSTNAME must be set")
        hec_handler = HECHandler(
            endpoint=os.environ["HEC_ENDPOINT"],
            token=os.environ["HEC_TOKEN"],
            hostname=hostname,
            source=os.environ["HEC_SOURCE"],
            sourcetype=os.environ["HEC_SOURCETYPE"],
        )
        hec_handler.setFormatter(formatter)
        logger.addHandler(hec_handler)

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
