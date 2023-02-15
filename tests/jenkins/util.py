from pathlib import Path


def create_paths(paths: list[str]) -> None:
    for path_str in paths:
        p = Path(path_str)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
