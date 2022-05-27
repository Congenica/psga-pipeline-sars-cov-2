#!/usr/bin/env python

from pathlib import Path

import click
from click import ClickException

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun
from scripts.db.queries import get_analysis_run
from scripts.util.logging import get_structlog_logger
from scripts.util.metadata import generate_notifications, inspect_metadata_file, ProcessedSamples

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)

NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS = {
    "illumina_artic": {"fastq", "bam"},
    "medaka_artic": {"fastq"},
    "no_ncov": {"fasta"},
}


def load_analysis_run_metadata(
    session,
    analysis_run_name: str,
    primer_scheme_name: str,
    primer_scheme_version: str,
    input_file_type: str,
    ncov_workflow: str,
    pipeline_version: str,
) -> AnalysisRun:
    """
    Set up a SARS-CoV-2 analysis run
    """
    analysis_run = get_analysis_run(session, analysis_run_name)

    if not analysis_run:
        analysis_run = AnalysisRun(
            analysis_run_name=analysis_run_name,
        )
        session.add(analysis_run)
    analysis_run.primer_scheme_name = primer_scheme_name
    analysis_run.primer_scheme_version = primer_scheme_version
    analysis_run.input_file_type = input_file_type.upper()
    analysis_run.ncov_workflow = ncov_workflow.upper()
    analysis_run.pipeline_version = pipeline_version

    return analysis_run


def validate(
    metadata_path: Path,
    analysis_run_name: str,
    primer_scheme_name: str,
    primer_scheme_version: str,
    input_file_type: str,
    ncov_workflow: str,
    pipeline_version: str,
    load_missing_samples: bool,
) -> ProcessedSamples:
    """
    Validate metadata and input parameters
    """

    if input_file_type not in NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS[ncov_workflow]:
        raise ClickException(f"ncov workflow '{ncov_workflow}' does not support input file type '{input_file_type}'")

    samples_with_two_reads = bool(ncov_workflow == "illumina_artic" and input_file_type == "fastq")

    with session_handler() as session:
        analysis_run = load_analysis_run_metadata(
            session,
            analysis_run_name,
            primer_scheme_name,
            primer_scheme_version,
            input_file_type,
            ncov_workflow,
            pipeline_version,
        )

        return inspect_metadata_file(session, metadata_path, analysis_run, samples_with_two_reads, load_missing_samples)


@click.command()
@click.option(
    "--metadata-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, readable=True),
    help="The metadata CSV input file",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
@click.option("--primer-scheme-name", required=True, type=str, help="The primer scheme name")
@click.option("--primer-scheme-version", required=True, type=str, help="The primer scheme version")
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
@click.option("--pipeline-version", type=str, required=True, help="mapping pipeline version")
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
@click.option(
    "--load-missing-samples",
    is_flag=True,
    help="This flag is only used for testing the pipeline in isolation."
    "In a non-testing execution, an error is raised if a sample is missing",
)
def check_metadata(
    metadata_path,
    analysis_run_name,
    primer_scheme_name,
    primer_scheme_version,
    input_file_type,
    ncov_workflow,
    pipeline_version,
    samples_with_valid_metadata_file,
    samples_with_invalid_metadata_file,
    load_missing_samples,
):
    """
    Read in a CSV file of metadata and check that samples are present in the database.
    If a sample in metadata is not found in the database, an error is raised,
    unless the flag --load-missing-samples is used.
    """

    samples = validate(
        Path(metadata_path),
        analysis_run_name,
        primer_scheme_name,
        primer_scheme_version,
        input_file_type,
        ncov_workflow,
        pipeline_version,
        load_missing_samples,
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
