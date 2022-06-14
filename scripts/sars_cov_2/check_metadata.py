#!/usr/bin/env python

from pathlib import Path

import click
from click import ClickException

from scripts.util.logging import get_structlog_logger
from scripts.util.metadata import generate_notifications, inspect_metadata_file, ProcessedSamples

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)

NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS = {
    "illumina_artic": {"fastq", "bam"},
    "medaka_artic": {"fastq"},
    "no_ncov": {"fasta"},
}


def validate(
    metadata_path: Path,
    input_file_type: str,
    ncov_workflow: str,
) -> ProcessedSamples:
    """
    Validate metadata and input parameters
    """

    if input_file_type not in NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS[ncov_workflow]:
        raise ClickException(f"ncov workflow '{ncov_workflow}' does not support input file type '{input_file_type}'")

    samples_with_two_reads = bool(ncov_workflow == "illumina_artic" and input_file_type == "fastq")

    return inspect_metadata_file(metadata_path, samples_with_two_reads)


@click.command()
@click.option(
    "--metadata-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, readable=True),
    help="The metadata CSV input file",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
@click.option(
    "--input-file-type",
    required=True,
    type=click.Choice(
        {ft for file_types in NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS.values() for ft in file_types},
        case_sensitive=True,
    ),
    help="The type of input files",
)
@click.option(
    "--ncov-workflow",
    required=True,
    type=click.Choice(set(NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS), case_sensitive=True),
    help="The name of the ncov workflow",
)
@click.option(
    "--samples-with-valid-metadata-file",
    required=True,
    type=click.Path(file_okay=True, writable=True),
    help="Output file path which will contain the samples with valid metadata",
)
@click.option(
    "--samples-with-invalid-metadata-file",
    required=True,
    type=click.Path(file_okay=True, writable=True),
    help="Output file path which will contain the samples with invalid metadata",
)
def check_metadata(
    metadata_path,
    analysis_run_name,
    input_file_type,
    ncov_workflow,
    samples_with_valid_metadata_file,
    samples_with_invalid_metadata_file,
):
    """
    Read a CSV metadata and check that is sound
    """

    samples = validate(
        Path(metadata_path),
        input_file_type,
        ncov_workflow,
    )

    generate_notifications(
        analysis_run_name,
        samples.valid,
        Path(samples_with_valid_metadata_file),
        samples.invalid,
        Path(samples_with_invalid_metadata_file),
    )

    if samples.invalid:
        raise ClickException("Errors encountered for sample ids: " + ", ".join(map(str, sorted(samples.invalid))))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    check_metadata()
