#!/usr/bin/env python

from pathlib import Path
import csv
import click
from click import ClickException

from scripts.util.logging import get_structlog_logger
from scripts.util.metadata import (
    ILLUMINA,
    ONT,
    UNKNOWN,
    SAMPLE_ID,
    SEQ_FILE_1,
    SEQ_FILE_2,
    EXPECTED_HEADERS,
    generate_notifications,
    normalise_row,
    ProcessedSamples,
)
from scripts.validation.check_csv_columns import check_csv_columns

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)


SEQUENCING_TECHNOLOGIES = [ILLUMINA, ONT, UNKNOWN]
FILE_NUM = "file_num"
SUPPORTED_FILES = {
    ILLUMINA: {
        "fastq": {FILE_NUM: 2},
        "fastq.gz": {FILE_NUM: 2},
        "fq": {FILE_NUM: 2},
        "fq.gz": {FILE_NUM: 2},
        "bam": {FILE_NUM: 1},
    },
    ONT: {
        "fastq": {FILE_NUM: 1},
        "fastq.gz": {FILE_NUM: 1},
        "fq": {FILE_NUM: 1},
        "fq.gz": {FILE_NUM: 1},
        "bam": {FILE_NUM: 1},
    },
    UNKNOWN: {
        "fasta": {FILE_NUM: 1},
        "fa": {FILE_NUM: 1},
    },
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
        supported_extensions = list(SUPPORTED_FILES[sequencing_technology])

        for row in reader:
            errs = []
            row = normalise_row(row)
            sample_id = row[SAMPLE_ID]

            # check sample_id
            if not sample_id:
                errs.append(f"{SAMPLE_ID} not available")

            # check file_1
            if not row[SEQ_FILE_1]:
                errs.append(f"{SEQ_FILE_1} for {sample_id} not available")
            else:
                file_1 = row[SEQ_FILE_1]
                # this returns 0 or 1 extension
                extensions = [ext for ext in supported_extensions if file_1.endswith(ext)]
                if not extensions:
                    errs.append(
                        f"{SAMPLE_ID}: {sample_id} has invalid file for sequencing technology {sequencing_technology}. "
                        f"Supported files are {supported_extensions}"
                    )
                else:
                    # check file_2 if required
                    file_1_ext = extensions[0]
                    two_reads = SUPPORTED_FILES[sequencing_technology][file_1_ext][FILE_NUM] == 2
                    if two_reads:
                        if not row[SEQ_FILE_2]:
                            errs.append(f"{SEQ_FILE_2} for {sample_id} not available")
                        elif not row[SEQ_FILE_2].endswith(file_1_ext):
                            errs.append(f"{SEQ_FILE_1} and {SEQ_FILE_2} for {sample_id} have different file types")

            if errs:
                sample_errors = "\n".join(errs)
                click.echo(f"Invalid row for {SAMPLE_ID} {sample_id}:\n{sample_errors}", err=True)
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
    type=click.Choice(SEQUENCING_TECHNOLOGIES, case_sensitive=True),
    help="the sequencer technology used for sequencing the samples",
)
def check_metadata(
    metadata_path,
    analysis_run_name,
    sequencing_technology,
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
        samples.invalid,
    )

    if samples.invalid:
        raise ClickException("Errors encountered for sample ids: " + ", ".join(map(str, sorted(samples.invalid))))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    check_metadata()
