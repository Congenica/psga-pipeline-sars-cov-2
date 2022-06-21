from pathlib import Path
import click

from scripts.util.logging import get_structlog_logger
from jenkins.sars_cov_2 import sars_cov_2

log_file = f"{Path(__file__).stem}.log"
get_structlog_logger(log_file=log_file)


@click.group()
@click.option(
    "--result-path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The calculated result file",
)
@click.option(
    "--expected-result-path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="The expected result file",
)
@click.option(
    "--psga-output-path",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    envvar="PSGA_OUTPUT_PATH",
    required=True,
    help="The PSGA pipeline path where all output files are stored",
)
@click.pass_context
def validate(ctx, result_path, expected_result_path, psga_output_path):
    """
    Compare the calculated result file against the expected result file.
    """
    ctx.ensure_object(dict)

    ctx.obj["result_path"] = Path(result_path)
    ctx.obj["expected_result_path"] = Path(expected_result_path)
    ctx.obj["psga_output_path"] = Path(psga_output_path)


validate.add_command(sars_cov_2)

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter,unexpected-keyword-arg
    validate(obj={})
