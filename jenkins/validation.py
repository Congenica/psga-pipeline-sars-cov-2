from pathlib import Path
import click

from scripts.util.logger import get_structlog_logger
from jenkins.sars_cov_2 import sars_cov_2

log_file = f"{Path(__file__).stem}.log"
get_structlog_logger(log_file=log_file)


@click.group()
@click.option(
    "--results-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The calculated result file",
)
@click.option(
    "--expected-results-csv",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The expected result file",
)
@click.option(
    "--output-path",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    required=True,
    help="The PSGA pipeline path where all output files are stored",
)
@click.pass_context
def validate(ctx, results_csv, expected_results_csv, output_path):
    """
    Compare the calculated result file against the expected result file.
    """
    ctx.ensure_object(dict)

    ctx.obj["results_csv"] = Path(results_csv)
    ctx.obj["expected_results_csv"] = Path(expected_results_csv)
    ctx.obj["output_path"] = output_path  # leave this as a string as it can be FS path or S3 uri


validate.add_command(sars_cov_2)

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter,unexpected-keyword-arg
    validate(obj={})
