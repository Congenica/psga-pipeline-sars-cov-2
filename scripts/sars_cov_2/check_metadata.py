#!/usr/bin/env python

from pathlib import Path
import csv
import click
from click import ClickException

from scripts.util.logging import get_structlog_logger
from scripts.util.metadata import (
    EXPECTED_HEADERS,
    generate_notifications,
    is_valid_uuid,
    normalise_row,
    ProcessedSamples,
)
from scripts.validation.check_csv_columns import check_csv_columns

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)

# the boolean set to True means 2 files (=reads), False means 1 file
SUPPORTED_FILES_BY_SEQUENCING_TECHNOLOGY = {
    "illumina": {"fastq.gz": True, "bam": False},
    "ont": {"fastq": False},
    "unknown": {"fasta": False},
}


def validate_metadata(
    metadata_path: Path,
    sequencing_technology: str,
) -> ProcessedSamples:
    """
    Validate the CSV metadata file of the samples to process
    """
    processed_samples = ProcessedSamples(
        valid=[],
        invalid=[],
    )

    with open(metadata_path, "r") as metadata_fd:
        reader = csv.DictReader(metadata_fd, delimiter=",")

        check_csv_columns(set(reader.fieldnames), EXPECTED_HEADERS)
        supported_extensions = list(SUPPORTED_FILES_BY_SEQUENCING_TECHNOLOGY[sequencing_technology].keys())

        for row in reader:
            errs = []
            row = normalise_row(row)
            sample_id = row["sample_id"]

            # check sample_id
            if not sample_id:
                errs.append("sample_id not available")
            elif not is_valid_uuid(sample_id):
                errs.append(f'sample_id "{sample_id}" is not a UUID')

            # check file_1
            if not row["file_1"]:
                errs.append(f"file_1 for {sample_id} not available")
            else:
                file_1 = row["file_1"]
                # this returns 0 or 1 extension
                extensions = [ext for ext in supported_extensions if file_1.endswith(ext)]
                if not extensions:
                    errs.append(
                        f"Sample: {sample_id} has invalid file for sequencing technology {sequencing_technology}. "
                        f"Supported files are {supported_extensions}"
                    )
                else:
                    # check file_2 if required
                    file_1_ext = extensions[0]
                    two_reads = (
                        sequencing_technology == "illumina"
                        and SUPPORTED_FILES_BY_SEQUENCING_TECHNOLOGY[sequencing_technology][file_1_ext]
                    )
                    if two_reads:
                        if not row["file_2"]:
                            errs.append(f"file_2 for {sample_id} not available")
                        elif not row["file_2"].endswith(file_1_ext):
                            errs.append(f"file_1 and file_2 for {sample_id} have different file types")

            if errs:
                sample_errors = "\n".join(errs)
                click.echo(f"Invalid row for sample {sample_id}:\n{sample_errors}", err=True)
                processed_samples.invalid.append(sample_id)
                continue

            processed_samples.valid.append(sample_id)

    return processed_samples


@click.command()
@click.option(
    "--metadata-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, readable=True),
    help="The metadata CSV input file",
)
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
@click.option(
    "--sequencing-technology",
    required=True,
    type=click.Choice(set(SUPPORTED_FILES_BY_SEQUENCING_TECHNOLOGY), case_sensitive=True),
    help="The name of the sequencing technology",
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
    sequencing_technology,
    samples_with_valid_metadata_file,
    samples_with_invalid_metadata_file,
):
    """
    Read a CSV metadata and check that is sound
    """

    samples = validate_metadata(
        Path(metadata_path),
        sequencing_technology,
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
