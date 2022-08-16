from pathlib import Path
from typing import List
import click

from jenkins.loading import load_data_from_csv, get_file_paths
from jenkins.compare import compare_merged_output_file, compare_output_files_set
from jenkins.config import data_config
from scripts.sars_cov_2.check_metadata import SEQUENCING_TECHNOLOGIES
from scripts.sars_cov_2.generate_pipeline_results_files import get_expected_output_files_per_sample, SampleIdResultFiles
from scripts.util.logging import get_structlog_logger
from scripts.util.metadata import UNKNOWN

logger = get_structlog_logger()


def get_expected_output_files(output_path: Path, sample_ids: List[str], sequencing_technology: str) -> List[Path]:
    """
    Return a list of ALL expected output paths
    """

    # In integration tests, we expect all samples to pass ncov qc if this is executed
    sample_ids_result_files = SampleIdResultFiles(
        all_samples=sample_ids,
        ncov_completed_samples=[] if sequencing_technology == UNKNOWN else sample_ids,
        ncov_qc_passed_samples=[] if sequencing_technology == UNKNOWN else sample_ids,
    )

    output_files_per_sample = get_expected_output_files_per_sample(
        output_path=output_path,
        sample_ids_result_files=sample_ids_result_files,
        sequencing_technology=sequencing_technology,
    )
    # generate a unified list of paths as non-sample results files must also be included
    output_files = [path for sample_paths in output_files_per_sample.values() for path in sample_paths]

    output_files.extend(
        [output_path / "logs" / f for f in ["check_metadata.log", "generate_pipeline_results_files.log"]]
    )

    notification_files = [
        "samples_failed_pangolin.txt",
        "samples_passed_pangolin.txt",
        "samples_unknown_pangolin.txt",
        "samples_with_invalid_metadata.txt",
        "samples_with_valid_metadata.txt",
    ]

    if sequencing_technology != UNKNOWN:
        notification_files.extend(
            [
                "samples_failed_ncov_qc.txt",
                "samples_passed_ncov_qc.txt",
                "samples_unknown_ncov_qc.txt",
            ]
        )
        output_files.append(output_path / "ncov2019-artic" / "ncov_qc.csv")

    output_files.append(output_path / "pangolin" / "all_lineages_report.csv")
    output_files.extend([output_path / "notifications" / p for p in notification_files])
    output_files.append(output_path / "results.csv")
    output_files.append(output_path / "resultfiles.json")

    return output_files


@click.command(name="sars_cov_2")
@click.option(
    "--sequencing-technology",
    required=True,
    type=click.Choice(SEQUENCING_TECHNOLOGIES, case_sensitive=True),
    help="The name of the sequencing technology",
)
@click.pass_context
def sars_cov_2(ctx, sequencing_technology: str):
    """
    Validate sars_cov_2 output
    """
    results_csv = ctx.obj["results_csv"]
    expected_results_csv = ctx.obj["expected_results_csv"]
    output_path = ctx.obj["output_path"]

    pathogen = "sars_cov_2"
    if pathogen not in data_config:
        ValueError(f"Configuration error. Pathogen {pathogen} not supported")
    data = data_config[pathogen]["config"]

    sample_ids = compare_merged_output_file(load_data_from_csv, data, results_csv, expected_results_csv)

    logger.info("Validation of output files set STARTED")
    exp_output_files = get_expected_output_files(output_path, sample_ids, sequencing_technology)
    calc_output_files = get_file_paths(output_path)
    compare_output_files_set(set(calc_output_files), set(exp_output_files))
    logger.info("Validation PASSED")
