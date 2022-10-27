from typing import Set
import csv
import click

from scripts.sars_cov_2.generate_pipeline_results_files import (
    EXPECTED_PRIMER_AUTODETECTION_HEADERS,
    EXPECTED_NCOV_HEADERS,
    EXPECTED_PANGOLIN_HEADERS,
)


def write_csv(filename: str, fieldnames: Set[str]) -> None:
    with open(filename, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()


@click.command()
def generate_default_files() -> None:
    """
    Generate empty CSV files for primer-autodetection, ncov and pangolin.
    These are used as default results, if the corresponding nextflow channel is empty.
    This script is used by the dockerfile. It is not used by the pipeline directly.

    This scripts was implemented so that the expected headers are not replicated
    """
    write_csv("primer_autodetection_empty.csv", EXPECTED_PRIMER_AUTODETECTION_HEADERS)
    write_csv("ncov_qc_empty.csv", EXPECTED_NCOV_HEADERS)
    write_csv("pangolin_empty.csv", EXPECTED_PANGOLIN_HEADERS)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_default_files()
