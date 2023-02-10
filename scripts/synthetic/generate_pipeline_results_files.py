from pathlib import Path
from typing import List, Set
from functools import partial, reduce
from itertools import cycle

import click
import pandas as pd

from scripts.util.logger import get_structlog_logger
from scripts.util.metadata import EXPECTED_HEADERS as EXPECTED_METADATA_HEADERS, SAMPLE_ID
from scripts.validation.check_csv_columns import check_csv_columns
from scripts.util.convert import csv_to_json
from scripts.util.data_loading import write_json
from scripts.util.slugs import get_file_with_type, FileType, RESULTFILES_TYPE

log_file = f"{Path(__file__).stem}.log"
logger = get_structlog_logger(log_file=log_file)


# These are pipeline generic columns
STATUS = "STATUS"

# Synthetic data from TB (fig 2 in https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00817-3)
PASS = "pass"
FAIL = "fail"
QC = "qc_status"
LINEAGE = "lineage"
RD = "region_of_difference"
SYNTHETIC_DATA = [
    {
        QC: PASS,
        LINEAGE: "1.1",
        RD: "RD239",
    },
    {
        QC: PASS,
        LINEAGE: "4.1",
        RD: "RD182",
    },
    {
        QC: PASS,
        LINEAGE: "2.1",
        RD: "RD105-RD207",
    },
    {
        QC: PASS,
        LINEAGE: "1.2",
        RD: "RD239",
    },
    {
        QC: PASS,
        LINEAGE: "2.2",
        RD: "RD150",
    },
    {
        QC: FAIL,
        LINEAGE: "3",
        RD: "RD750",
    },
    {
        QC: FAIL,
        LINEAGE: "4.1",
        RD: "RD193",
    },
    {
        QC: PASS,
        LINEAGE: "5",
        RD: "RD711",
    },
]


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


def _generate_results_csv(
    all_samples: List[str],
    df_synthetic: pd.DataFrame,
    qc_unrelated_failing_samples: List[str],
    output_results_csv_file: str,
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
    dfs_to_merge = [df_status, df_synthetic]

    # partial stores part of a function’s arguments resulting in a new object with a simplified signature.
    # reduce applies cumulatively the new partial object to the items of iterable (list of dataframes here).
    merge = partial(pd.merge, how="outer")
    df_merged = reduce(merge, dfs_to_merge)

    # upper case column header
    df_merged.columns = [col.upper() for col in df_merged.columns]
    # move sample_id col to first column
    sample_id_col_data = df_merged.pop(SAMPLE_ID)
    df_merged.insert(0, SAMPLE_ID, sample_id_col_data)

    # save to CSV
    df_merged.to_csv(output_results_csv_file, encoding="utf-8", index=False)


def get_expected_output_files_per_sample(
    output_path: str,
    sample_ids_result_files: List[str],
) -> RESULTFILES_TYPE:
    """
    Return a dictionary {sample_id, list_of_expected_output_paths}
    """
    # initialise the dictionary keys
    output_files: RESULTFILES_TYPE = {sample_id: [] for sample_id in sample_ids_result_files}

    for sample_id in sample_ids_result_files:
        output_files[sample_id] = get_file_with_type(
            output_path=output_path,
            inner_dirs=["analysis"],
            filetypes=[FileType("_analysis.txt", "txt/analysis")],
            sample_id=sample_id,
        )

    return output_files


def _generate_resultfiles_json(
    sample_ids_result_files: List[str],
    output_path: str,
    output_resultfiles_json_file: Path,
) -> None:
    """
    Generate a JSON file containing the list of expected result files per sample
    """
    output_files_per_sample = get_expected_output_files_per_sample(
        output_path=output_path,
        sample_ids_result_files=sample_ids_result_files,
    )

    write_json(output_files_per_sample, output_resultfiles_json_file)


@click.command()
@click.option(
    "--metadata-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="the sample metadata file",
)
@click.option(
    "--output-results-csv-file",
    type=click.Path(dir_okay=False, writable=True),
    default="results.csv",
    help="CSV file containing the summary of results",
)
@click.option(
    "--output-results-json-file",
    type=click.Path(dir_okay=False, writable=True),
    default="results.json",
    help="JSON file containing the summary of results. This is results.csv but in JSON format.",
)
@click.option(
    "--output-resultfiles-json-file",
    type=click.Path(dir_okay=False, writable=True),
    default="resultfiles.json",
    help="JSON file containing all the expected files per sample",
)
@click.option(
    "--output-path",
    type=str,
    required=True,
    help="output_path output path where sample result files are stored (e.g. s3://bucket/path/analysis_run)",
)
def generate_pipeline_results_files(
    metadata_file: str,
    output_results_csv_file: str,
    output_results_json_file: str,
    output_resultfiles_json_file: str,
    output_path: str,
) -> None:
    """
    Generate pipeline results files
    """
    # generate a cycle list so that generated data for each sample is predictable
    SYNTHETIC_DATA_POOL = cycle(SYNTHETIC_DATA)

    df_metadata = load_data_from_csv(Path(metadata_file), EXPECTED_METADATA_HEADERS)
    all_samples = sorted(df_metadata[SAMPLE_ID].tolist())

    # generate some synthetic results.
    successful_samples = all_samples
    qc_unrelated_failing_samples: List[str] = []

    synthetic_sample_results = [next(SYNTHETIC_DATA_POOL) for sample in successful_samples]

    synthetic_results = {
        SAMPLE_ID: successful_samples,
        QC: [s[QC] for s in synthetic_sample_results],
        LINEAGE: [s[LINEAGE] for s in synthetic_sample_results],
        RD: [s[RD] for s in synthetic_sample_results],
    }
    df_synthetic = pd.DataFrame(data=synthetic_results)

    _generate_results_csv(all_samples, df_synthetic, qc_unrelated_failing_samples, output_results_csv_file)

    csv_to_json(output_results_csv_file, output_results_json_file, SAMPLE_ID)

    _generate_resultfiles_json(successful_samples, output_path, Path(output_resultfiles_json_file))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_pipeline_results_files()
