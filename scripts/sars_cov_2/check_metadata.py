#!/usr/bin/env python

from typing import List
import csv
from datetime import datetime
from pathlib import Path

import click
from click import ClickException

from sqlalchemy.orm import scoped_session
from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample
from scripts.util.metadata import METADATA_FILE_EXPECTED_HEADERS, generate_notifications, is_valid_uuid


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
    analysis_run = (
        session.query(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )
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


def validate_and_normalise_row(row, ncov_workflow, input_file_type):

    # strip leading and trailing spaces from everything
    for f in METADATA_FILE_EXPECTED_HEADERS:
        row[f] = row[f].lstrip().rstrip() if row[f] is not None else ""

    errs = []

    # add validations here
    sample_id = row["sample_id"]
    if not sample_id:
        errs.append("sample_id not available")
    elif not is_valid_uuid(sample_id):
        errs.append(f'sample_id "{sample_id}" is not a UUID')

    if not row["file_1"]:
        errs.append(f"file_1 for {sample_id} not available")
    if not row["md5_1"]:
        errs.append(f"md5_1 for {sample_id} not available")

    if ncov_workflow == "illumina_artic" and input_file_type == "fastq":
        if not row["file_2"]:
            errs.append(f"file_2 for {sample_id} not available")
        if not row["md5_2"]:
            errs.append(f"md5_2 for {sample_id} not available")

    if errs:
        raise ValueError("\n".join(errs))

    return row


def validate(
    session: scoped_session,
    metadata_path: str,
    analysis_run_name: str,
    primer_scheme_name: str,
    primer_scheme_version: str,
    input_file_type: str,
    ncov_workflow: str,
    pipeline_version: str,
    load_missing_samples: bool,
) -> tuple[List[str], List[str]]:
    """
    Validate metadata and input parameters
    """

    if input_file_type not in NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS[ncov_workflow]:
        raise ClickException(f"ncov workflow '{ncov_workflow}' does not support input file type '{input_file_type}'")

    reader = csv.DictReader(metadata_path, delimiter=",")

    headers = set(reader.fieldnames)
    if not METADATA_FILE_EXPECTED_HEADERS.issubset(headers):
        err = (
            "Unexpected CSV headers, got:\n"
            + ", ".join(headers)
            + "\n, but expect at least \n"
            + ", ".join(METADATA_FILE_EXPECTED_HEADERS)
        )
        raise ClickException(err)

    valid_samples = []
    invalid_samples = []
    # only if --load-missing-samples is used
    loaded_samples = []

    analysis_run = load_analysis_run_metadata(
        session,
        analysis_run_name,
        primer_scheme_name,
        primer_scheme_version,
        input_file_type,
        ncov_workflow,
        pipeline_version,
    )

    for row in reader:
        try:
            row = validate_and_normalise_row(row, ncov_workflow, input_file_type)
            sample_name = row["sample_id"]

            existing_sample = (
                session.query(Sample)
                .join(AnalysisRun)
                .filter(
                    Sample.sample_name == sample_name,
                    AnalysisRun.analysis_run_name == analysis_run_name,
                )
                .one_or_none()
            )

            if existing_sample:
                existing_sample.analysis_run_id = analysis_run.analysis_run_id
            elif load_missing_samples:
                # the pipeline metadata does not contain the date_collected field
                yesterday = datetime.strftime(datetime.now(), "%Y-%m-%d")
                sample = Sample(
                    analysis_run_id=analysis_run.analysis_run_id,
                    sample_name=sample_name,
                    date_collected=yesterday,
                    metadata_loaded=True,
                )
                session.add(sample)
                loaded_samples.append(sample_name)
            else:
                raise ValueError(f"Sample {sample_name} not found in the database, but listed in pipeline metadata")

            session.commit()
            valid_samples.append(sample_name)

        except ValueError as e:
            sample_name = row["sample_id"]
            click.echo(f"Invalid row for sample {sample_name}:\n{e}", err=True)
            invalid_samples.append(sample_name)

    if load_missing_samples:
        click.echo("Inserted samples: " + ", ".join(map(str, loaded_samples)))

    return valid_samples, invalid_samples


@click.command()
@click.option("--metadata-path", required=True, type=click.File("r"), help="The metadata CSV input file")
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

    with session_handler() as session:
        valid_samples, invalid_samples = validate(
            session,
            metadata_path,
            analysis_run_name,
            primer_scheme_name,
            primer_scheme_version,
            input_file_type,
            ncov_workflow,
            pipeline_version,
            load_missing_samples,
        )

    generate_notifications(
        valid_samples,
        Path(samples_with_valid_metadata_file),
        invalid_samples,
        Path(samples_with_invalid_metadata_file),
    )

    if invalid_samples:
        raise ClickException("Errors encountered for sample ids: " + ", ".join(map(str, sorted(invalid_samples))))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    check_metadata()
