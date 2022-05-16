#!/usr/bin/env python

from typing import List
import csv
from datetime import datetime
from pathlib import Path
from uuid import UUID

import click
from click import ClickException

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample
from scripts.util.data_dumping import write_list_to_file

METADATA_FILE_EXPECTED_HEADERS = {"sample_id", "file_1", "file_2", "md5_1", "md5_2"}
INPUT_FILE_TYPES = {"bam", "fastq"}
WORKFLOWS = {"illumina_artic", "medaka_artic"}


class MedakaArticWithBamError(Exception):
    pass


def is_valid_uuid(uuid_str: str) -> bool:
    try:
        UUID(uuid_str)
        return True
    except ValueError:
        return False


def load_analysis_run_metadata(
    session,
    analysis_run_name: str,
    primer_scheme_name: str,
    primer_scheme_version: str,
    input_file_type: str,
    workflow: str,
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
    analysis_run.workflow = workflow.upper()
    analysis_run.pipeline_version = pipeline_version

    return analysis_run


def _validate_and_normalise_row(row, workflow, input_file_type):

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

    if workflow == "illumina_artic" and input_file_type == "fastq":
        if not row["file_2"]:
            errs.append(f"file_2 for {sample_id} not available")
        if not row["md5_2"]:
            errs.append(f"md5_2 for {sample_id} not available")

    if errs:
        raise ValueError("\n".join(errs))

    return row


def write_sample_list_files(
    valid_samples: List[str],
    file_current_samples: str,
) -> None:
    if file_current_samples:
        write_list_to_file(valid_samples, Path(file_current_samples))


@click.command()
@click.option("--file", required=True, type=click.File("r"), help="The metadata CSV input file")
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
@click.option("--primer-scheme-name", required=True, type=str, help="The primer scheme name")
@click.option("--primer-scheme-version", required=True, type=str, help="The primer scheme version")
@click.option(
    "--input-file-type",
    required=True,
    type=click.Choice(INPUT_FILE_TYPES, case_sensitive=True),
    help="The type of input files",
)
@click.option(
    "--workflow", required=True, type=click.Choice(WORKFLOWS, case_sensitive=True), help="The name of the workflow"
)
@click.option("--pipeline-version", type=str, required=True, help="mapping pipeline version")
@click.option(
    "--output-current-samples-with-metadata",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with sample names within the current metadata file",
)
@click.option(
    "--load-missing-samples",
    is_flag=True,
    help="This flag is only used for testing the pipeline in isolation."
    "In a non-testing execution, an error is raised if a sample is missing",
)
def check_metadata(
    file,
    analysis_run_name,
    output_current_samples_with_metadata,
    primer_scheme_name,
    primer_scheme_version,
    input_file_type,
    workflow,
    pipeline_version,
    load_missing_samples,
):
    """
    Read in a CSV file of metadata and check that samples are present in the database.
    If a sample in metadata is not found in the database, an error is raised,
    unless the flag --load-missing-samples is used.
    """

    if input_file_type == "bam" and workflow == "medaka_artic":
        click.echo("Error: medaka_artic workflow does not support input bam files")
        raise MedakaArticWithBamError

    reader = csv.DictReader(file, delimiter=",")

    headers = set(reader.fieldnames)
    if not METADATA_FILE_EXPECTED_HEADERS.issubset(headers):
        err = (
            "Unexpected CSV headers, got:\n"
            + ", ".join(headers)
            + "\n, but expect at least \n"
            + ", ".join(METADATA_FILE_EXPECTED_HEADERS)
        )
        raise ClickException(err)

    file_samples_with_valid_meta = []

    # only if --load-missing-samples is used
    inserted = set()

    errors = set()
    with session_handler() as session:

        analysis_run = load_analysis_run_metadata(
            session,
            analysis_run_name,
            primer_scheme_name,
            primer_scheme_version,
            input_file_type,
            workflow,
            pipeline_version,
        )

        for row in reader:
            try:
                row = _validate_and_normalise_row(row, workflow, input_file_type)
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
                    inserted.add(sample_name)
                else:
                    raise ValueError(f"Sample {sample_name} not found in the database, but listed in pipeline metadata")

                session.commit()
                file_samples_with_valid_meta.append(sample_name)

            except ValueError as e:
                click.echo(f"Invalid row for sample ID {row['sample_id']}:\n{e}", err=True)
                errors.add(row["sample_id"])

    write_sample_list_files(
        valid_samples=file_samples_with_valid_meta,
        file_current_samples=output_current_samples_with_metadata,
    )

    if errors:
        raise ClickException("Errors encountered for sample ids: " + ", ".join(map(str, sorted(errors))))

    if load_missing_samples:
        click.echo("Inserted samples: " + ", ".join(map(str, inserted)))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    check_metadata()
