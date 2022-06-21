from pathlib import Path
from typing import List
import click

from jenkins.loading import load_data_from_csv, get_file_paths
from jenkins.compare import compare_merged_output_file, compare_output_files_set
from jenkins.config import data_config
from scripts.util.logging import get_structlog_logger
from scripts.sars_cov_2.check_metadata import NCOV_WORKFLOW_FILE_TYPE_VALID_COMBINATIONS

logger = get_structlog_logger()


def get_expected_output_files(root: Path, sample_names: List[str], filetype: str, ncov_workflow: str) -> List[Path]:
    """
    Return a list of expected output paths
    """
    output_files: List[Path] = []

    output_files.extend([root / "logs" / f for f in ["check_metadata.log", "merge_ncov_pangolin_csv_files.log"]])
    output_files.append(root / "merged_output" / "pipeline_output.csv")
    output_files.append(root / "pangolin" / "all_lineages_report.csv")
    output_files.extend([root / "reheadered-fasta" / f"{s}.fasta" for s in sample_names])

    notification_files = [
        "samples_failed_pangolin.txt",
        "samples_passed_pangolin.txt",
        "samples_unknown_pangolin.txt",
        "samples_with_invalid_metadata.txt",
        "samples_with_valid_metadata.txt",
    ]

    if ncov_workflow != "no_ncov":
        notification_files.extend(
            ["samples_failed_ncov_qc.txt", "samples_passed_ncov_qc.txt", "samples_unknown_ncov_qc.txt"]
        )

        if ncov_workflow == "illumina_artic" and filetype == "fastq":
            fastqc_suffixes = [f"{r}_fastqc.{e}" for r in [1, 2] for e in ["html", "zip"]]
        else:
            fastqc_suffixes = [f"fastqc.{e}" for e in ["html", "zip"]]
        output_files.extend([root / "fastqc" / f"{s}_{e}" for s in sample_names for e in fastqc_suffixes])

        if ncov_workflow == "illumina_artic":
            ncov_fasta_suffixes = [".primertrimmed.consensus.fa"]
            ncov_plots_suffixes = [".depth.png"]
        elif ncov_workflow == "medaka_artic":
            ncov_fasta_suffixes = [".consensus.fa", ".muscle.in.fa", ".muscle.out.fa", ".preconsensus.fa"]
            ncov_plots_suffixes = ["-barplot.png", "-boxplot.png", ".depth.png"]
        output_files.extend(
            [root / "ncov2019-artic" / "output_fasta" / f"{s}{e}" for s in sample_names for e in ncov_fasta_suffixes]
        )
        output_files.extend(
            [root / "ncov2019-artic" / "output_plots" / f"{s}{e}" for s in sample_names for e in ncov_plots_suffixes]
        )
        output_files.append(root / "ncov2019-artic" / "ncov_qc.csv")

    output_files.extend([root / "notifications" / p for p in notification_files])

    return output_files


@click.command(name="sars_cov_2")
@click.option(
    "--filetype",
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
@click.pass_context
def sars_cov_2(ctx, filetype: str, ncov_workflow: str):
    """
    Validate sars_cov_2 output
    """
    result_path = ctx.obj["result_path"]
    expected_result_path = ctx.obj["expected_result_path"]
    psga_output_path = ctx.obj["psga_output_path"]

    pathogen = "sars_cov_2"
    if pathogen not in data_config:
        ValueError(f"Configuration error. Pathogen {pathogen} not supported")
    data = data_config[pathogen]["config"]

    sample_names = compare_merged_output_file(load_data_from_csv, data, result_path, expected_result_path)

    logger.info("Validation of output files set STARTED")
    exp_output_files = get_expected_output_files(psga_output_path, sample_names, filetype, ncov_workflow)
    calc_output_files = get_file_paths(psga_output_path)
    compare_output_files_set(set(calc_output_files), set(exp_output_files))
    logger.info("Validation PASSED")
