from dataclasses import dataclass, field
from marshmallow import validate

from app.scripts.util.logger import LOG_LEVELS, get_structlog_logger


@dataclass
class Event:
    analysis_run: str = field(metadata={"required": True})
    level: str = field(metadata={"required": True, "validate": validate.OneOf(LOG_LEVELS)})
    message: str = field(metadata={"required": True})
    samples: list[str] = field(metadata={"required": True}, default_factory=list)


@dataclass
class Notification:
    events: dict[str, Event] = field(metadata={"required": True}, default_factory=dict)

    def publish(self) -> None:
        """
        Publish all the notifications to file and as log messages
        """
        for evt in self.events.values():
            samples = evt.samples

            # log per sample
            logger_level = getattr(get_structlog_logger(), evt.level.lower())
            for sample in samples:
                logger_level(evt.message, sample=sample, analysis_run=evt.analysis_run)
