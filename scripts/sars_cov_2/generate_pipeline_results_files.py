from pathlib import Path, PosixPath
from typing import Dict, List, Set, Tuple
from dataclasses import dataclass, field
from functools import partial, reduce
import json
from json import JSONEncoder

import click
import pandas as pd

from scripts.util.logging import get_structlog_logger, ERROR, WARNING, INFO
from scripts.util.metadata import EXPECTED_HEADERS as EXPECTED_METADATA_HEADERS, SAMPLE_ID, ILLUMINA, ONT, UNKNOWN
from scripts.util.notifications import Event, Notification
from scripts.validation.check_csv_columns import check_csv_columns

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)


# header of ncov qc summary CSV file
NCOV_SAMPLE_ID_COL = "sample_name"
EXPECTED_NCOV_HEADERS = {
    NCOV_SAMPLE_ID_COL,
    "pct_N_bases",
    "pct_covered_bases",
    "longest_no_N_run",
    "num_aligned_reads",
    "qc_pass",
    "fasta",
    "bam",
}

# header of pangolin lineages CSV file
PANGOLIN_SAMPLE_ID_COL = "taxon"
EXPECTED_PANGOLIN_HEADERS = {
    PANGOLIN_SAMPLE_ID_COL,
    "lineage",
    "conflict",
    "ambiguity_score",
    "scorpio_call",
    "scorpio_support",
    "scorpio_conflict",
    "scorpio_notes",
    "version",
    "pangolin_version",
    "scorpio_version",
    "constellation_version",
    "is_designated",
    "qc_status",
    "qc_notes",
    "note",
}


# These are pipeline generic columns
STATUS = "STATUS"

# these columns point to specific files and are not needed
COLUMNS_TO_REMOVE_FROM_RESULTS_CSV = {
    "FASTA",
    "BAM",
}

# notification event names
UNKNOWN_NCOV = "unknown_ncov"
FAILED_NCOV = "failed_ncov"
PASSED_NCOV = "passed_ncov"
UNKNOWN_PANGOLIN = "unknown_pangolin"
FAILED_PANGOLIN = "failed_pangolin"
PASSED_PANGOLIN = "passed_pangolin"

# output files storing sample ids for unknown/failing/passed ncov/pangolin QC
SAMPLES_UNKNOWN_NCOV_QC_FILE = f"samples_{UNKNOWN_NCOV}_qc.txt"
SAMPLES_FAILED_NCOV_QC_FILE = f"samples_{FAILED_NCOV}_qc.txt"
SAMPLES_PASSED_NCOV_QC_FILE = f"samples_{PASSED_NCOV}_qc.txt"
SAMPLES_UNKNOWN_PANGOLIN_FILE = f"samples_{UNKNOWN_PANGOLIN}.txt"
SAMPLES_FAILED_PANGOLIN_FILE = f"samples_{FAILED_PANGOLIN}.txt"
SAMPLES_PASSED_PANGOLIN_FILE = f"samples_{PASSED_PANGOLIN}.txt"


@dataclass
class SampleIdResultFiles:
    """
    Organisation of the expected results files per sample
    """

    # all samples of the analysis run
    all_samples: List[str] = field(metadata={"required": True}, default_factory=list)
    # samples which completed ncov, whether passing or failing ncov qc
    ncov_completed_samples: List[str] = field(metadata={"required": True}, default_factory=list)
    # samples passing ncov qc
    ncov_qc_passed_samples: List[str] = field(metadata={"required": True}, default_factory=list)


class PathJSONEncoder(JSONEncoder):
    """
    Enable JSON serialisation of PosixPath objects
    """

    def default(self, o):
        if isinstance(o, PosixPath):
            return str(o)
        return super().default(o)


def load_data_from_csv(
    csv_path: Path, expected_columns: Set[str], sample_name_col_to_rename: str = None
) -> pd.DataFrame:
    """
    Load the CSV content to a Pandas dataframe. An arbitrary column name used to indentify the sample id
    can be renamed to "sample_id"
    """
    df = pd.read_csv(csv_path)
    check_csv_columns(set(df.columns), expected_columns)
    if sample_name_col_to_rename:
        df = df.rename(columns={sample_name_col_to_rename: SAMPLE_ID})
    return df


def _generate_notifications(
    analysis_run_name: str,
    all_samples: List[str],
    df_ncov: pd.DataFrame,
    df_pangolin: pd.DataFrame,
    notifications_path: Path,
) -> Tuple[List[str], Notification]:
    """
    Generate and publish output pipeline notifications.
    Return the list of samples which failed not due to QC
    """
    events = {}

    # initialise ncov as if it had not executed. This is the default case in which fastas were processed
    ncov_all_samples = all_samples
    qc_unrelated_failing_ncov_samples = []
    ncov_samples_passing_qc = all_samples

    if not df_ncov.empty:
        ncov_all_samples = df_ncov[SAMPLE_ID].tolist()
        qc_unrelated_failing_ncov_samples = [s for s in all_samples if s not in ncov_all_samples]
        ncov_samples_failing_qc = df_ncov.loc[~df_ncov["qc_pass"]][SAMPLE_ID]
        ncov_samples_passing_qc = df_ncov.loc[df_ncov["qc_pass"]][SAMPLE_ID]

        events = {
            UNKNOWN_NCOV: Event(
                analysis_run=analysis_run_name,
                path=Path(notifications_path / SAMPLES_UNKNOWN_NCOV_QC_FILE),
                level=ERROR,
                message="ncov QC unknown",
                samples=qc_unrelated_failing_ncov_samples,
            ),
            FAILED_NCOV: Event(
                analysis_run=analysis_run_name,
                path=Path(notifications_path / SAMPLES_FAILED_NCOV_QC_FILE),
                level=WARNING,
                message="ncov QC failed",
                samples=ncov_samples_failing_qc,
            ),
            PASSED_NCOV: Event(
                analysis_run=analysis_run_name,
                path=Path(notifications_path / SAMPLES_PASSED_NCOV_QC_FILE),
                level=INFO,
                message="ncov QC passed",
                samples=ncov_samples_passing_qc,
            ),
        }

    # pangolin is executed after ncov, unless the latter is skipped (e.g. fasta files)
    pangolin_all_samples = df_pangolin[SAMPLE_ID].tolist()
    # NOTE: pangolin is not executed for samples failing ncov QC
    qc_unrelated_failing_pangolin_samples = [s for s in ncov_samples_passing_qc if s not in pangolin_all_samples]
    pangolin_samples_failing_qc = df_pangolin.loc[df_pangolin["qc_status"] == "fail"][SAMPLE_ID]
    pangolin_samples_passing_qc = df_pangolin.loc[df_pangolin["qc_status"] == "pass"][SAMPLE_ID]

    events = {
        **events,
        UNKNOWN_PANGOLIN: Event(
            analysis_run=analysis_run_name,
            path=Path(notifications_path / SAMPLES_UNKNOWN_PANGOLIN_FILE),
            level=ERROR,
            message="pangolin QC unknown",
            samples=qc_unrelated_failing_pangolin_samples,
        ),
        FAILED_PANGOLIN: Event(
            analysis_run=analysis_run_name,
            path=Path(notifications_path / SAMPLES_FAILED_PANGOLIN_FILE),
            level=WARNING,
            message="pangolin QC failed",
            samples=pangolin_samples_failing_qc,
        ),
        PASSED_PANGOLIN: Event(
            analysis_run=analysis_run_name,
            path=Path(notifications_path / SAMPLES_PASSED_PANGOLIN_FILE),
            level=INFO,
            message="pangolin QC passed",
            samples=pangolin_samples_passing_qc,
        ),
    }

    notifications = Notification(events=events)
    notifications.publish()

    qc_unrelated_failing_samples = list(
        set(qc_unrelated_failing_ncov_samples) | set(qc_unrelated_failing_pangolin_samples)
    )

    return qc_unrelated_failing_samples, events


def _generate_results_csv(
    all_samples: List[str],
    df_ncov: pd.DataFrame,
    df_pangolin: pd.DataFrame,
    qc_unrelated_failing_samples: List[str],
    output_csv_file: str,
) -> None:
    """
    Generate the pipeline results CSV file.
    """

    status_col_data = [
        {SAMPLE_ID: s, STATUS: "Failed"} if s in qc_unrelated_failing_samples else {SAMPLE_ID: s, STATUS: "Completed"}
        for s in all_samples
    ]
    df_status = pd.DataFrame(status_col_data)

    # list of dataframes to merge. They share 1 shared column: SAMPLE_ID
    dfs_to_merge = [df_status, df_ncov, df_pangolin]

    # partial stores part of a function’s arguments resulting in a new object with a simplified signature.
    # reduce applies cumulatively the new partial object to the items of iterable (list of dataframes here).
    merge = partial(pd.merge, how="outer")
    df_merged = reduce(merge, dfs_to_merge)

    # upper case column header
    df_merged.columns = [col.upper() for col in df_merged.columns]
    # move sample_id col to first column
    sample_id_col_data = df_merged.pop(SAMPLE_ID)
    df_merged.insert(0, SAMPLE_ID, sample_id_col_data)

    # remove unwanted columns
    df_merged.drop(COLUMNS_TO_REMOVE_FROM_RESULTS_CSV, axis=1, inplace=True)
    # save to CSV
    df_merged.to_csv(output_csv_file, encoding="utf-8", index=False)


def get_expected_output_files_per_sample(
    output_path: Path,
    sample_ids_result_files: SampleIdResultFiles,
    sequencing_technology: str,
) -> Dict[str, List[Path]]:
    """
    Return a dictionary {sample_id, list_of_expected_output_paths}
    """
    # initialise the dictionary keys
    output_files: Dict[str, List[Path]] = {sample_id: [] for sample_id in sample_ids_result_files.all_samples}

    def _append_reheadered_fasta(output_path, sample_id, output_files):
        output_files[sample_id].append(output_path / "reheadered-fasta" / f"{sample_id}.fasta")

    if sequencing_technology == UNKNOWN:
        # FASTA samples
        for sample_id in sample_ids_result_files.all_samples:
            # expected files for all samples
            _append_reheadered_fasta(output_path, sample_id, output_files)
    else:
        if sequencing_technology == ILLUMINA:
            contamination_removal_suffixes = [f"_{r}.fastq.gz" for r in [1, 2]]
            fastqc_suffixes = [f"{r}_fastqc.zip" for r in [1, 2]]
            ncov_bam_suffixes = [".mapped.primertrimmed.sorted.bam", ".mapped.primertrimmed.sorted.bam.bai"]
            ncov_fasta_suffixes = [".primertrimmed.consensus.fa"]
            ncov_variants_suffixes = [".variants.tsv"]
        elif sequencing_technology == ONT:
            contamination_removal_suffixes = [".fastq.gz"]
            fastqc_suffixes = ["fastqc.zip"]
            ncov_bam_suffixes = [".primertrimmed.rg.sorted.bam", ".primertrimmed.rg.sorted.bam.bai"]
            ncov_fasta_suffixes = [".consensus.fa", ".muscle.in.fa", ".muscle.out.fa", ".preconsensus.fa"]
            ncov_variants_suffixes = [".pass.vcf.gz", ".pass.vcf.gz.tbi"]
        else:
            raise ValueError(f"Unsupported sequencing_technology: {sequencing_technology}")

        for sample_id in sample_ids_result_files.all_samples:
            # expected files for all samples
            output_files[sample_id].extend(
                [
                    output_path / "contamination_removal" / "cleaned_fastq" / f"{sample_id}{e}"
                    for e in contamination_removal_suffixes
                ]
            )
            output_files[sample_id].append(output_path / "contamination_removal" / "counting" / f"{sample_id}.txt")
            output_files[sample_id].extend([output_path / "fastqc" / f"{sample_id}_{e}" for e in fastqc_suffixes])

        for sample_id in sample_ids_result_files.ncov_completed_samples:
            # expected files for samples which completed ncov
            output_files[sample_id].extend(
                [output_path / "ncov2019-artic" / "output_bam" / f"{sample_id}{e}" for e in ncov_bam_suffixes]
            )
            output_files[sample_id].extend(
                [output_path / "ncov2019-artic" / "output_fasta" / f"{sample_id}{e}" for e in ncov_fasta_suffixes]
            )
            output_files[sample_id].extend(
                [output_path / "ncov2019-artic" / "output_variants" / f"{sample_id}{e}" for e in ncov_variants_suffixes]
            )
            output_files[sample_id].append(output_path / "ncov2019-artic" / "output_plots" / f"{sample_id}.depth.png")

        for sample_id in sample_ids_result_files.ncov_qc_passed_samples:
            # expected files for samples which passed ncov QC
            _append_reheadered_fasta(output_path, sample_id, output_files)

    return output_files


def _generate_resultfiles_json(
    sequencing_technology: str,
    sample_ids_result_files: SampleIdResultFiles,
    output_path: Path,
    output_json_file: Path,
) -> None:
    """
    Generate a JSON file containing the list of expected result files per sample
    """
    output_files_per_sample = get_expected_output_files_per_sample(
        output_path=output_path,
        sample_ids_result_files=sample_ids_result_files,
        sequencing_technology=sequencing_technology,
    )

    with open(output_json_file, "w") as outfile:
        json.dump(output_files_per_sample, outfile, cls=PathJSONEncoder, sort_keys=True, indent=4)


@click.command()
@click.option("--analysis-run-name", required=True, type=str, help="The name of the analysis run")
@click.option(
    "--metadata-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="the sample metadata file",
)
@click.option(
    "--ncov-qc-csv-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    help="ncov pipeline resulting qc csv file",
)
@click.option(
    "--pangolin-csv-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="pangolin pipeline resulting csv file",
)
@click.option(
    "--output-csv-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file merging the results from ncov and pangolin",
)
@click.option(
    "--output-json-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file containing all the expected files per sample",
)
@click.option(
    "--output-path",
    type=str,
    required=True,
    help="output_path output path where sample result files are stored (e.g. s3://bucket/path/analysis_run)",
)
@click.option(
    "--sequencing-technology",
    type=click.Choice([ILLUMINA, ONT, UNKNOWN], case_sensitive=True),
    required=True,
    help="the sequencer technology used for sequencing the samples",
)
@click.option(
    "--notifications-path",
    type=click.Path(dir_okay=True, writable=True),
    default=".",
    help="path used for saving the output notification files. This is a directory",
)
def generate_pipeline_results_files(
    analysis_run_name: str,
    metadata_file: str,
    ncov_qc_csv_file: str,
    pangolin_csv_file: str,
    output_csv_file: str,
    output_json_file: str,
    output_path: str,
    sequencing_technology: str,
    notifications_path: str,
) -> None:
    """
    Generate pipeline results files
    """
    df_metadata = load_data_from_csv(Path(metadata_file), EXPECTED_METADATA_HEADERS)
    all_samples = df_metadata[SAMPLE_ID].tolist()

    if ncov_qc_csv_file:
        # ncov was executed
        df_ncov = load_data_from_csv(Path(ncov_qc_csv_file), EXPECTED_NCOV_HEADERS, NCOV_SAMPLE_ID_COL)
    else:
        # create a dataframe with header but no rows
        df_ncov = pd.DataFrame({c: [] for c in EXPECTED_NCOV_HEADERS})
        df_ncov = df_ncov.rename(columns={NCOV_SAMPLE_ID_COL: SAMPLE_ID})

    df_pangolin = load_data_from_csv(Path(pangolin_csv_file), EXPECTED_PANGOLIN_HEADERS, PANGOLIN_SAMPLE_ID_COL)

    qc_unrelated_failing_samples, events = _generate_notifications(
        analysis_run_name, all_samples, df_ncov, df_pangolin, Path(notifications_path)
    )

    _generate_results_csv(all_samples, df_ncov, df_pangolin, qc_unrelated_failing_samples, output_csv_file)

    sample_ids_result_files = SampleIdResultFiles(
        all_samples=all_samples,
        ncov_completed_samples=[]
        if sequencing_technology == UNKNOWN or not {FAILED_NCOV, PASSED_NCOV} <= events.keys()
        else list(set(events[FAILED_NCOV].samples) | set(events[PASSED_NCOV].samples)),
        ncov_qc_passed_samples=[]
        if sequencing_technology == UNKNOWN or PASSED_NCOV not in events
        else list(events[PASSED_NCOV].samples),
    )

    _generate_resultfiles_json(
        sequencing_technology, sample_ids_result_files, Path(output_path), Path(output_json_file)
    )


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_pipeline_results_files()
