from pathlib import Path
from typing import List


def create_paths(paths: List[str]) -> None:
    for path_str in paths:
        p = Path(path_str)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
