from pathlib import Path
from typing import Dict, List
from dataclasses import dataclass, field
from marshmallow import validate

from scripts.util.data_dumping import write_list_to_file
from scripts.util.logging import LOG_LEVELS, get_structlog_logger


@dataclass
class Event:
    analysis_run: str = field(metadata={"required": True})
    path: Path = field(metadata={"required": True})
    level: str = field(metadata={"required": True, "validate": validate.OneOf(LOG_LEVELS)})
    message: str = field(metadata={"required": True})
    samples: List[str] = field(metadata={"required": True}, default_factory=list)


@dataclass
class Notification:
    events: Dict[str, Event] = field(metadata={"required": True}, default_factory=dict)

    def publish(self):
        """
        Publish all the notifications to file and as log messages
        """
        for evt in self.events.values():
            samples = evt.samples

            # log per sample as in the requirements
            logger_level = getattr(get_structlog_logger(), evt.level.lower())
            for sample in samples:
                logger_level(evt.message, sample=sample, analysis_run=evt.analysis_run)

            # write samples to file
            write_list_to_file(samples, evt.path)
