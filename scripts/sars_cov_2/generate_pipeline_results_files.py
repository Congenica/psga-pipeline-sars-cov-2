from pathlib import Path
from typing import List, Set
from functools import partial, reduce
import click
import pandas as pd

from scripts.util.logging import get_structlog_logger
from scripts.util.metadata import EXPECTED_HEADERS as EXPECTED_METADATA_HEADERS
from scripts.util.notifications import Event, Notification
from scripts.validation.check_csv_columns import check_csv_columns

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)

EXPECTED_NCOV_HEADERS = {
    "sample_name",
    "pct_N_bases",
    "pct_covered_bases",
    "longest_no_N_run",
    "num_aligned_reads",
    "qc_pass",
    "fasta",
    "bam",
}
EXPECTED_PANGOLIN_HEADERS = {
    "taxon",
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
SAMPLE_ID_COL = "SAMPLE_ID"
STATUS = "STATUS"

# these columns point to specific files and are not needed
COLUMNS_TO_REMOVE_FROM_RESULTS_CSV = {
    "FASTA",
    "BAM",
}

# output files storing sample ids for unknown/failing/passed ncov/pangolin QC
SAMPLES_UNKNOWN_NCOV_QC_FILE = "samples_unknown_ncov_qc.txt"
SAMPLES_FAILED_NCOV_QC_FILE = "samples_failed_ncov_qc.txt"
SAMPLES_PASSED_NCOV_QC_FILE = "samples_passed_ncov_qc.txt"
SAMPLES_UNKNOWN_PANGOLIN_FILE = "samples_unknown_pangolin.txt"
SAMPLES_FAILED_PANGOLIN_FILE = "samples_failed_pangolin.txt"
SAMPLES_PASSED_PANGOLIN_FILE = "samples_passed_pangolin.txt"


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
        df = df.rename(columns={sample_name_col_to_rename: SAMPLE_ID_COL})
    return df


def _generate_notifications(
    analysis_run_name: str,
    all_samples: List[str],
    df_ncov: pd.DataFrame,
    df_pangolin: pd.DataFrame,
    notifications_path: Path,
) -> List[str]:
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
        ncov_all_samples = df_ncov[SAMPLE_ID_COL].tolist()
        qc_unrelated_failing_ncov_samples = [s for s in all_samples if s not in ncov_all_samples]
        ncov_samples_failing_qc = df_ncov.loc[~df_ncov["qc_pass"]][SAMPLE_ID_COL]
        ncov_samples_passing_qc = df_ncov.loc[df_ncov["qc_pass"]][SAMPLE_ID_COL]

        events = {
            "unknown_ncov": Event(
                analysis_run=analysis_run_name,
                path=Path(notifications_path / SAMPLES_UNKNOWN_NCOV_QC_FILE),
                level="ERROR",
                message="ncov QC unknown",
                samples=qc_unrelated_failing_ncov_samples,
            ),
            "failed_ncov": Event(
                analysis_run=analysis_run_name,
                path=Path(notifications_path / SAMPLES_FAILED_NCOV_QC_FILE),
                level="WARNING",
                message="ncov QC failed",
                samples=ncov_samples_failing_qc,
            ),
            "passed_ncov": Event(
                analysis_run=analysis_run_name,
                path=Path(notifications_path / SAMPLES_PASSED_NCOV_QC_FILE),
                level="INFO",
                message="ncov QC passed",
                samples=ncov_samples_passing_qc,
            ),
        }

    # pangolin is executed after ncov, unless the latter is skipped (e.g. fasta files)
    pangolin_all_samples = df_pangolin[SAMPLE_ID_COL].tolist()
    # NOTE: pangolin is not executed for samples failing ncov QC
    qc_unrelated_failing_pangolin_samples = [s for s in ncov_samples_passing_qc if s not in pangolin_all_samples]
    pangolin_samples_failing_qc = df_pangolin.loc[df_pangolin["qc_status"] == "fail"][SAMPLE_ID_COL]
    pangolin_samples_passing_qc = df_pangolin.loc[df_pangolin["qc_status"] == "pass"][SAMPLE_ID_COL]

    events = {
        **events,
        "unknown_pangolin": Event(
            analysis_run=analysis_run_name,
            path=Path(notifications_path / SAMPLES_UNKNOWN_PANGOLIN_FILE),
            level="ERROR",
            message="pangolin QC unknown",
            samples=qc_unrelated_failing_pangolin_samples,
        ),
        "failed_pangolin": Event(
            analysis_run=analysis_run_name,
            path=Path(notifications_path / SAMPLES_FAILED_PANGOLIN_FILE),
            level="WARNING",
            message="pangolin QC failed",
            samples=pangolin_samples_failing_qc,
        ),
        "passed_pangolin": Event(
            analysis_run=analysis_run_name,
            path=Path(notifications_path / SAMPLES_PASSED_PANGOLIN_FILE),
            level="INFO",
            message="pangolin QC passed",
            samples=pangolin_samples_passing_qc,
        ),
    }

    notifications = Notification(events=events)
    notifications.publish()

    qc_unrelated_failing_samples = list(
        set(qc_unrelated_failing_ncov_samples) | set(qc_unrelated_failing_pangolin_samples)
    )

    return qc_unrelated_failing_samples


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
        {SAMPLE_ID_COL: s, STATUS: "Failed"}
        if s in qc_unrelated_failing_samples
        else {SAMPLE_ID_COL: s, STATUS: "Completed"}
        for s in all_samples
    ]
    df_status = pd.DataFrame(status_col_data)

    # list of dataframes to merge. They share 1 shared column: SAMPLE_ID_COL
    dfs_to_merge = [df_status, df_ncov, df_pangolin]

    # partial stores part of a function’s arguments resulting in a new object with a simplified signature.
    # reduce applies cumulatively the new partial object to the items of iterable (list of dataframes here).
    merge = partial(pd.merge, how="outer")
    df_merged = reduce(merge, dfs_to_merge)

    # upper case column header
    df_merged.columns = [col.upper() for col in df_merged.columns]
    # move sample_id col to first column
    sample_id_col_data = df_merged.pop(SAMPLE_ID_COL)
    df_merged.insert(0, SAMPLE_ID_COL, sample_id_col_data)

    # remove unwanted columns
    df_merged.drop(COLUMNS_TO_REMOVE_FROM_RESULTS_CSV, axis=1, inplace=True)
    # save to CSV
    df_merged.to_csv(output_csv_file, encoding="utf-8", index=False)


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
    notifications_path: str,
) -> None:
    """
    Generate pipeline results files
    """
    df_metadata = load_data_from_csv(Path(metadata_file), EXPECTED_METADATA_HEADERS)
    all_samples = df_metadata[SAMPLE_ID_COL].tolist()

    if ncov_qc_csv_file:
        # ncov was executed
        df_ncov = load_data_from_csv(Path(ncov_qc_csv_file), EXPECTED_NCOV_HEADERS, "sample_name")
    else:
        # create a dataframe with header but no rows
        df_ncov = pd.DataFrame({c: [] for c in EXPECTED_NCOV_HEADERS})
        df_ncov = df_ncov.rename(columns={"sample_name": SAMPLE_ID_COL})

    df_pangolin = load_data_from_csv(Path(pangolin_csv_file), EXPECTED_PANGOLIN_HEADERS, "taxon")

    # generate notifications
    qc_unrelated_failing_samples = _generate_notifications(
        analysis_run_name, all_samples, df_ncov, df_pangolin, Path(notifications_path)
    )

    _generate_results_csv(all_samples, df_ncov, df_pangolin, qc_unrelated_failing_samples, output_csv_file)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_pipeline_results_files()
