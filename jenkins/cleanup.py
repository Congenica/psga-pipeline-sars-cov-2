from pathlib import Path
import click

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample
from scripts.util.logging import get_structlog_logger

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)


def remove_analysis_run(analysis_run_name: str) -> None:
    """
    Remove an analysis run and all the associated samples from the database.
    The DB schema and model are not set up to delete on cascade for safety reasons,
    therefore here we delete data explicitly
    """
    with session_handler() as session:
        samples = (
            session.query(Sample)
            .join(AnalysisRun)
            .filter(
                AnalysisRun.analysis_run_name == analysis_run_name,
            )
            .all()
        )

        if not samples:
            logger.info(f"Analysis run {analysis_run_name} not found in the DB. Nothing to delete")
            return

        sample_names = [s.sample_name for s in samples]
        logger.info(f"Deleting analysis run-associated samples: {sample_names}")
        for s in samples:
            session.delete(s)

        logger.info(f"Deleting analysis_run: {analysis_run_name}")
        session.query(AnalysisRun).filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        ).delete()


@click.command()
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
def cleanup(analysis_run_name):
    """
    Clean up the database so that it is read for running a new integration test.
    """
    logger.info(f"Cleaning up analysis run: {analysis_run_name}")
    remove_analysis_run(analysis_run_name)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    cleanup()
