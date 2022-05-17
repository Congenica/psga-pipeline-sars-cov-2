from dataclasses import dataclass, field
from typing import Dict, List
from pathlib import Path
from marshmallow import validate

from scripts.util.data_dumping import write_list_to_file


LOG_LEVELS = {"ERROR", "WARNING", "INFO", "DEBUG"}


@dataclass
class Event:
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
            # log this event

            # write to file
            write_list_to_file(evt["samples"], evt["path"])
