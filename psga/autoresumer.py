#!/usr/bin/env python3

import time
from pathlib import Path
import subprocess
from typing import Dict
import uuid
import click
import psutil

SESSION_ID_DIR_ENV_VAR = "PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH"
MAX_ATTEMPTS_ENV_VAR = "PSGA_MAX_ATTEMPTS"
SLEEP_TIME_BETWEEN_ATTEMPTS_ENV_VAR = "PSGA_SLEEP_TIME_BETWEEN_ATTEMPTS"  # seconds

SESSION_ID_TYPE = Dict[str, str]


def is_uuid(val: str) -> bool:
    try:
        uuid.UUID(val)
        return True
    except ValueError:
        return False


def get_session_ids(session_id_dir: Path) -> SESSION_ID_TYPE:
    """
    Return a dictionary of session_id, run_id.
    <session_id> is the Nextflow session-id
    """
    # look for files like: <run_id>_<session_id>, where <session_id> is a uuid
    session_ids: SESSION_ID_TYPE = {}
    files = [f.stem for f in Path(session_id_dir).iterdir() if f.is_file() and f.stem.count("_") > 0 and not f.suffix]

    for f in files:
        # session_id does not contain underscores, but run_id could.
        # here we split reverse and stop at the first '_' found
        run_id, session_id = f.rsplit("_", 1)
        if is_uuid(session_id):
            session_ids[session_id] = run_id

    return session_ids


def is_process_running(process_name: str, program_name: str) -> bool:
    """
    Check if there is any running process that contains the given name process_name.
    """
    # Iterate over the all the running process
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if process_name in proc.name().lower() and program_name in proc.cmdline():
                return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False


def resume_nextflow(session_id_dir: Path, max_attempts: int) -> None:
    """
    Resume a nextflow session if present
    """

    session_ids = get_session_ids(session_id_dir)

    for session_id, run_id in session_ids.items():
        stem = f"{run_id}_{session_id}"
        attempt = max(
            [int(f.suffix[1:]) for f in Path(session_id_dir).iterdir() if f.is_file() and f.stem == stem and f.suffix],
            default=0,
        )

        if attempt < max_attempts:
            new_attempt = attempt + 1
            print(f"Resuming Nextflow pipeline: run_id: {run_id}, session_id: {session_id} - attempt #: {new_attempt}")
            # note: <run_id>_<session_id>.0 does not exist
            Path(session_id_dir / f"{stem}.{attempt}").unlink(missing_ok=True)
            Path(session_id_dir / f"{stem}.{new_attempt}").touch()

            session_id_path = Path(session_id_dir / stem)
            # run the command. If the command fails, this script does not care intentionally
            subprocess.run(["/bin/bash", str(session_id_path), session_id], check=False)
            break


@click.command()
@click.option(
    "--session-id-dir",
    type=Path,
    envvar=SESSION_ID_DIR_ENV_VAR,
    default="/data",
    required=True,
    help="The directory containing the session id files (if any)",
)
@click.option(
    "--max-attempts",
    type=int,
    envvar=MAX_ATTEMPTS_ENV_VAR,
    default=3,
    required=True,
    help="The maximum number of attempts before discarding the pipeline resumation",
)
@click.option(
    "--sleep-time",
    type=int,
    envvar=SLEEP_TIME_BETWEEN_ATTEMPTS_ENV_VAR,
    default=60,
    required=True,
    help="The sleep time between attemtps in seconds",
)
def autoresumer(session_id_dir, max_attempts, sleep_time):
    """
    This script attempts to resume a Nextflow session for 3 times. It runs automatically in the background.
    """
    session_id_path = Path(session_id_dir)
    session_id_path.mkdir(parents=True, exist_ok=True)

    while True:
        if not is_process_running("java", "nextflow.cli.Launcher"):
            resume_nextflow(session_id_path, max_attempts)
        time.sleep(sleep_time)


if __name__ == "__main__":
    autoresumer()  # pylint: disable=no-value-for-parameter,unexpected-keyword-arg
