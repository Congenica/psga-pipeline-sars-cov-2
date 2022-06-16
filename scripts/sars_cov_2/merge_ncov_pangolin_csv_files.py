from pathlib import Path
from typing import Set
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

# these columns point to specific files and are not needed
COLUMNS_TO_REMOVE_FROM_MERGED_OUTPUT = {
    "fasta",
    "bam",
}


def load_data_from_csv(
    csv_path: Path, expected_columns: Set[str], sample_name_col_to_rename: str = "sample_id"
) -> pd.DataFrame:
    """
    Load the CSV content to a Pandas dataframe
    """
    df = pd.read_csv(csv_path)
    check_csv_columns(set(df.columns), expected_columns)
    df = df.rename(columns={sample_name_col_to_rename: "sample_id"})
    return df


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
    "--merged-output-csv-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file merging the results from ncov and pangolin",
)
@click.option(
    "--samples-unknown-ncov-qc-file",
    type=click.Path(dir_okay=False, writable=True),
    default="samples_unknown_ncov_qc.txt",
    help="output file storing the samples with unknown ncov QC",
)
@click.option(
    "--samples-failed-ncov-qc-file",
    type=click.Path(dir_okay=False, writable=True),
    default="samples_failed_ncov_qc.txt",
    help="output file storing the samples failing ncov QC",
)
@click.option(
    "--samples-passed-ncov-qc-file",
    type=click.Path(dir_okay=False, writable=True),
    default="samples_passed_ncov_qc.txt",
    help="output file storing the samples passing ncov QC",
)
@click.option(
    "--samples-unknown-pangolin-file",
    type=click.Path(dir_okay=False, writable=True),
    default="samples_unknown_pangolin.txt",
    help="output file storing the samples with unknown pangolin QC",
)
@click.option(
    "--samples-failed-pangolin-file",
    type=click.Path(dir_okay=False, writable=True),
    default="samples_failed_pangolin.txt",
    help="output file storing the samples failing pangolin QC",
)
@click.option(
    "--samples-passed-pangolin-file",
    type=click.Path(dir_okay=False, writable=True),
    default="samples_passed_pangolin.txt",
    help="output file storing the samples passing pangolin QC",
)
def merge_ncov_pangolin_csv_files(
    analysis_run_name: str,
    metadata_file: str,
    ncov_qc_csv_file: str,
    pangolin_csv_file: str,
    merged_output_csv_file: str,
    samples_unknown_ncov_qc_file: str,
    samples_failed_ncov_qc_file: str,
    samples_passed_ncov_qc_file: str,
    samples_unknown_pangolin_file: str,
    samples_failed_pangolin_file: str,
    samples_passed_pangolin_file: str,
) -> None:
    """
    Merge ncov and pangolin CSV files
    """
    sample_id = "sample_id"

    df_metadata = load_data_from_csv(Path(metadata_file), EXPECTED_METADATA_HEADERS, sample_id)

    if ncov_qc_csv_file:
        # ncov was executed
        df_ncov = load_data_from_csv(Path(ncov_qc_csv_file), EXPECTED_NCOV_HEADERS, "sample_name")
    else:
        # create a dataframe with header but no rows
        df_ncov = pd.DataFrame({c: [] for c in EXPECTED_NCOV_HEADERS})
        df_ncov = df_ncov.rename(columns={"sample_name": "sample_id"})

    df_pangolin = load_data_from_csv(Path(pangolin_csv_file), EXPECTED_PANGOLIN_HEADERS, "taxon")

    df_merged = pd.merge(df_ncov, df_pangolin, how="outer")

    # move sample_id col to first column
    sample_id_col = df_merged.pop(sample_id)
    df_merged.insert(0, sample_id, sample_id_col)

    # remove unwanted columns
    df_merged.drop(COLUMNS_TO_REMOVE_FROM_MERGED_OUTPUT, axis=1, inplace=True)

    df_merged.to_csv(merged_output_csv_file, encoding="utf-8", index=False)

    # generate notifications
    all_samples = df_metadata[sample_id].tolist()
    ncov_all_samples = df_ncov[sample_id].tolist()
    pangolin_all_samples = df_pangolin[sample_id].tolist()

    events = {}

    if not df_ncov.empty:
        events = {
            "unknown_ncov": Event(
                analysis_run=analysis_run_name,
                path=samples_unknown_ncov_qc_file,
                level="ERROR",
                message="ncov QC unknown",
                samples=[s for s in all_samples if s not in ncov_all_samples],
            ),
            "failed_ncov": Event(
                analysis_run=analysis_run_name,
                path=samples_failed_ncov_qc_file,
                level="WARNING",
                message="ncov QC failed",
                samples=df_ncov.loc[~df_ncov["qc_pass"]][sample_id],
            ),
            "passed_ncov": Event(
                analysis_run=analysis_run_name,
                path=samples_passed_ncov_qc_file,
                level="INFO",
                message="ncov QC passed",
                samples=df_ncov.loc[df_ncov["qc_pass"]][sample_id],
            ),
        }
    else:
        # analysis run is a batch of fasta sample, ncov did not run
        ncov_all_samples = all_samples

    # pangolin is executed after ncov, unless this is skipped

    events = {
        **events,
        "unknown_pangolin": Event(
            analysis_run=analysis_run_name,
            path=samples_unknown_pangolin_file,
            level="ERROR",
            message="pangolin QC unknown",
            samples=[s for s in ncov_all_samples if s not in pangolin_all_samples],
        ),
        "failed_pangolin": Event(
            analysis_run=analysis_run_name,
            path=samples_failed_pangolin_file,
            level="WARNING",
            message="pangolin QC failed",
            samples=df_pangolin.loc[df_pangolin["qc_status"] == "fail"][sample_id],
        ),
        "passed_pangolin": Event(
            analysis_run=analysis_run_name,
            path=samples_passed_pangolin_file,
            level="INFO",
            message="pangolin QC passed",
            samples=df_pangolin.loc[df_pangolin["qc_status"] == "pass"][sample_id],
        ),
    }

    notifications = Notification(events=events)
    notifications.publish()


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    merge_ncov_pangolin_csv_files()
