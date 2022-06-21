from pathlib import Path
from typing import List


def create_paths(paths: List[Path]) -> None:
    for p in paths:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
