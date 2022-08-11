from pathlib import Path
from typing import List
import click

from jenkins.loading import load_data_from_csv, get_file_paths
from jenkins.compare import compare_merged_output_file, compare_output_files_set
from jenkins.config import data_config
from scripts.util.logging import get_structlog_logger
from scripts.sars_cov_2.check_metadata import SEQUENCING_TECHNOLOGIES

logger = get_structlog_logger()


def get_expected_output_files(root: Path, sample_ids: List[str], sequencing_technology: str) -> List[Path]:
    """
    Return a list of expected output paths
    """
    output_files: List[Path] = []

    output_files.extend([root / "logs" / f for f in ["check_metadata.log", "merge_ncov_pangolin_csv_files.log"]])

    notification_files = [
        "samples_failed_pangolin.txt",
        "samples_passed_pangolin.txt",
        "samples_unknown_pangolin.txt",
        "samples_with_invalid_metadata.txt",
        "samples_with_valid_metadata.txt",
    ]

    if sequencing_technology != "unknown":
        # skip ncov
        notification_files.extend(
            ["samples_failed_ncov_qc.txt", "samples_passed_ncov_qc.txt", "samples_unknown_ncov_qc.txt"]
        )

        if sequencing_technology == "illumina":
            contamination_removal_suffixes = [f"_{r}.fastq.gz" for r in [1, 2]]
            fastqc_suffixes = [f"{r}_fastqc.zip" for r in [1, 2]]
            ncov_fasta_suffixes = [".primertrimmed.consensus.fa"]
        elif sequencing_technology == "ont":
            contamination_removal_suffixes = [".fastq.gz"]
            fastqc_suffixes = ["fastqc.zip"]
            ncov_fasta_suffixes = [".consensus.fa", ".muscle.in.fa", ".muscle.out.fa", ".preconsensus.fa"]
        else:
            raise ValueError(f"Unsupported sequencing_technology: {sequencing_technology}")

        for sample_id in sample_ids:
            output_files.extend(
                [
                    root / "contamination_removal" / "cleaned_fastq" / f"{sample_id}{e}"
                    for e in contamination_removal_suffixes
                ]
            )
            output_files.append(root / "contamination_removal" / "counting" / f"{sample_id}.txt")
            output_files.extend([root / "fastqc" / f"{sample_id}_{e}" for e in fastqc_suffixes])
            output_files.extend(
                [root / "ncov2019-artic" / "output_fasta" / f"{sample_id}{e}" for e in ncov_fasta_suffixes]
            )
            output_files.append(root / "ncov2019-artic" / "output_plots" / f"{sample_id}.depth.png")
            output_files.append(root / "reheadered-fasta" / f"{sample_id}.fasta")

        output_files.append(root / "ncov2019-artic" / "ncov_qc.csv")

    output_files.append(root / "pangolin" / "all_lineages_report.csv")
    output_files.append(root / "merged_output" / "pipeline_output.csv")
    output_files.extend([root / "notifications" / p for p in notification_files])

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
    result_path = ctx.obj["result_path"]
    expected_result_path = ctx.obj["expected_result_path"]
    psga_output_path = ctx.obj["psga_output_path"]

    pathogen = "sars_cov_2"
    if pathogen not in data_config:
        ValueError(f"Configuration error. Pathogen {pathogen} not supported")
    data = data_config[pathogen]["config"]

    sample_ids = compare_merged_output_file(load_data_from_csv, data, result_path, expected_result_path)

    logger.info("Validation of output files set STARTED")
    exp_output_files = get_expected_output_files(psga_output_path, sample_ids, sequencing_technology)
    calc_output_files = get_file_paths(psga_output_path)
    compare_output_files_set(set(calc_output_files), set(exp_output_files))
    logger.info("Validation PASSED")
