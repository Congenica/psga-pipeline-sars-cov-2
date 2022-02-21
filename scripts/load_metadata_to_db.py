#!/usr/bin/env python

from typing import List
import csv
from datetime import date
import re
from pathlib import Path

import click
from click import ClickException

from scripts.db.database import session_handler
from scripts.db.models import AnalysisRun, Sample, SampleQC

METADATA_FILE_EXPECTED_HEADERS = {"SAMPLE ID", "ASSIGN DATE"}
INPUT_FILE_TYPES = {"bam", "fastq"}
WORKFLOWS = {"illumina_artic", "medaka_artic"}


class MedakaArticWithBamError(Exception):
    pass


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
    analysis_run.input_file_type = input_file_type
    analysis_run.workflow = workflow
    analysis_run.pipeline_version = pipeline_version

    return analysis_run


def _validate_and_normalise_row(row):

    # strip leading and trailing spaces from everything
    for f in METADATA_FILE_EXPECTED_HEADERS:
        row[f] = row[f].lstrip().rstrip() if row[f] is not None else ""

    errs = []

    # ASSIGN DATE should be dd/mm/yyyy
    assign_date_match = re.match(r"(\d{2})/(\d{2})/(\d{4})$", row["ASSIGN DATE"])
    if assign_date_match:
        try:
            row["ASSIGN DATE"] = date(
                int(assign_date_match.group(3)), int(assign_date_match.group(2)), int(assign_date_match.group(1))
            )
        except ValueError as e:
            errs.append(f"ASSIGN DATE \"{row['ASSIGN DATE']}\" is not a valid date: {e}")
    else:
        errs.append(f"ASSIGN DATE \"{row['ASSIGN DATE']}\" doesn't look like a date")

    if errs:
        raise ValueError("\n".join(errs))

    return row


def write_list_to_file(out_file: Path, elements: List[str]) -> None:
    with open(out_file, "w") as f:
        for element in elements:
            f.write(f"{element}\n")


def write_sample_list_files(
    updated_samples: List[str],
    valid_samples: List[str],
    file_all_samples: str,
    file_current_samples: str,
    file_updated_samples: str,
    file_qc_pass_samples: str,
) -> None:
    if file_all_samples:
        with session_handler() as session:
            all_samples_with_metadata = session.query(Sample.sample_name).filter(Sample.metadata_loaded).all()
            all_samples_with_metadata = [x[0] for x in all_samples_with_metadata]
            write_list_to_file(Path(file_all_samples), all_samples_with_metadata)
    if file_current_samples:
        write_list_to_file(Path(file_current_samples), valid_samples)
    if file_qc_pass_samples:
        with session_handler() as session:
            samples_with_qc_pass = (
                session.query(Sample.sample_name).join(Sample.sample_qc).filter(SampleQC.qc_pass).all()
            )
            samples_with_qc_pass = [x[0] for x in samples_with_qc_pass]
            write_list_to_file(Path(file_qc_pass_samples), samples_with_qc_pass)
    if file_updated_samples:
        write_list_to_file(Path(file_updated_samples), updated_samples)


@click.command()
@click.option("--file", required=True, type=click.File("r"), help="The metadata TSV input file")
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
    "--output-all-samples-with-metadata",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with all ALL sample names found historically, "
    "which have metadata loaded in the database",
)
@click.option(
    "--output-current-samples-with-metadata",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with sample names within the current metadata file, "
    "which were loaded into the database",
)
@click.option(
    "--output-samples-with-qc-pass",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with all sample names, "
    "which were processed by artic-ncov2019 pipeline and were marked as QC_PASS",
)
@click.option(
    "--output-samples-updated",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with all sample names, which were re-submitted and overwritten with provided metadata",
)
def load_metadata(
    file,
    analysis_run_name,
    output_all_samples_with_metadata,
    output_current_samples_with_metadata,
    output_samples_with_qc_pass,
    output_samples_updated,
    primer_scheme_name,
    primer_scheme_version,
    input_file_type,
    workflow,
    pipeline_version,
):
    """
    Read in a TSV file of metadata and load each row into the database. Invalid rows are warned about, but
    skipped over.
    """

    if input_file_type == "bam" and workflow == "medaka_artic":
        click.echo("Error: medaka_artic workflow does not support input bam files")
        raise MedakaArticWithBamError

    reader = csv.DictReader(file, delimiter="\t")

    headers = set(reader.fieldnames)
    if not METADATA_FILE_EXPECTED_HEADERS.issubset(headers):
        err = (
            "Unexpected TSV headers, got:\n"
            + ", ".join(headers)
            + "\n, but expect at least \n"
            + ", ".join(METADATA_FILE_EXPECTED_HEADERS)
        )
        raise ClickException(err)

    file_samples_with_valid_meta = []
    file_samples_updated = []

    inserted = set()
    updated = set()
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
                row = _validate_and_normalise_row(row)
                sample_name = row["SAMPLE ID"]

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
                    existing_sample.date_collected = row["ASSIGN DATE"]
                    existing_sample.metadata_loaded = True
                    updated.add(sample_name)
                    file_samples_updated.append(sample_name)
                else:
                    sample = Sample(
                        analysis_run_id=analysis_run.analysis_run_id,
                        sample_name=sample_name,
                        date_collected=row["ASSIGN DATE"],
                        metadata_loaded=True,
                    )
                    session.add(sample)
                    inserted.add(sample_name)
                session.commit()
                file_samples_with_valid_meta.append(sample_name)
            except ValueError as e:
                click.echo(f"Invalid row for sample ID {row['SAMPLE ID']}:\n{e}", err=True)
                errors.add(row["SAMPLE ID"])

    write_sample_list_files(
        updated_samples=file_samples_updated,
        valid_samples=file_samples_with_valid_meta,
        file_all_samples=output_all_samples_with_metadata,
        file_current_samples=output_current_samples_with_metadata,
        file_updated_samples=output_samples_updated,
        file_qc_pass_samples=output_samples_with_qc_pass,
    )

    if errors:
        raise ClickException("Errors encountered: " + ", ".join(map(str, errors)))
    click.echo("Inserted samples: " + ", ".join(map(str, inserted)))
    click.echo("Updated samples: " + ", ".join(map(str, updated)))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_metadata()
