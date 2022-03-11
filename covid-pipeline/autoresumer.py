#!/usr/bin/env python3

import time
from pathlib import Path
import subprocess
import uuid
import click
import psutil

SESSION_ID_DIR_ENV_VAR = "COVID_PIPELINE_OUTPUT_PATH"
MAX_ATTEMPTS_ENV_VAR = "MAX_ATTEMPTS"
MAX_ATTEMPTS = 3
SLEEP_TIME = 60  # seconds


def is_uuid(val: str) -> bool:
    try:
        uuid.UUID(val)
        return True
    except ValueError:
        return False


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
    session_ids = [f.stem for f in Path(session_id_dir).iterdir() if f.is_file() and is_uuid(f.stem) and not f.suffix]

    for session_id in session_ids:
        attempt = max(
            [
                int(f.suffix[1:])
                for f in Path(session_id_dir).iterdir()
                if f.is_file() and f.stem == session_id and f.suffix
            ],
            default=0,
        )

        if attempt < max_attempts:
            new_attempt = attempt + 1
            print(f"Resuming Nextflow pipeline. Attempt: {new_attempt}. Session ID: {session_id}")
            # note: <session_id>.0 does not exist
            Path(session_id_dir / f"{session_id}.{attempt}").unlink(missing_ok=True)
            Path(session_id_dir / f"{session_id}.{new_attempt}").touch()

            session_id_path = Path(session_id_dir / session_id)
            # run the command. If the command fails, this script does not care intentionally
            subprocess.run(["/bin/bash", str(session_id_path), session_id], check=False)
            break


@click.command()
@click.option(
    "--session-id-dir",
    type=click.Path(exists=True, resolve_path=True),
    envvar=SESSION_ID_DIR_ENV_VAR,
    default="/data/input",
    required=True,
    help="The directory containing the session id files (if any)",
)
@click.option(
    "--max-attempts",
    type=int,
    envvar=MAX_ATTEMPTS_ENV_VAR,
    default=MAX_ATTEMPTS,
    required=True,
    help="The maximum number of attempts before discarding the pipeline resumation",
)
def autoresumer(session_id_dir, max_attempts):
    """
    This script attempts to resume a Nextflow session for 3 times. It runs automatically in the background.
    """
    session_id_path = Path(session_id_dir)
    while True:
        if not is_process_running("java", "nextflow.cli.Launcher"):
            resume_nextflow(session_id_path, max_attempts)
        time.sleep(SLEEP_TIME)


if __name__ == "__main__":
    autoresumer()  # pylint: disable=no-value-for-parameter,unexpected-keyword-arg
