from pathlib import Path
from typing import Tuple
import pandas as pd
import click

from concat_csv import concat

SAMTOOLS_COVERAGE_TSV_EXTENSION = ".coverage.tsv"
PRIMER_DETECTION_SUFFIX = "_primer_detection.csv"
PRIMER_DATA_SUFFIX = "_primer_data.csv"

RNAME_COL = "#rname"
PRIMER_AUTODETECTION_SAMPLE_ID_COL = "sample_id"
PRIMER_AUTODETECTION_PRIMER_INPUT_COL = "primer_input"
PRIMER_AUTODETECTION_PRIMER_COL = "primer_detected"
PRIMER_AUTODETECTION_NUMREADS_COL = "primer_numreads"
PRIMER_AUTODETECTION_COVBASES_COL = "primer_covbases"
PRIMER_AUTODETECTION_COVERAGE_COL = "primer_coverage"
EXPECTED_PRIMER_AUTODETECTION_HEADERS = {
    PRIMER_AUTODETECTION_SAMPLE_ID_COL,
    PRIMER_AUTODETECTION_PRIMER_INPUT_COL,
    PRIMER_AUTODETECTION_PRIMER_COL,
    "startpos",
    "endpos",
    PRIMER_AUTODETECTION_NUMREADS_COL,
    PRIMER_AUTODETECTION_COVBASES_COL,
    PRIMER_AUTODETECTION_COVERAGE_COL,
    "meandepth",
    "meanbaseq",
    "meanmapq",
}
PRIMER_AUTODETECTION_PRIMER_SCORE_COL = PRIMER_AUTODETECTION_NUMREADS_COL
# columns to be renamed
SAMTOOLS_COVERAGE_TO_PSGA_COLS = {
    RNAME_COL: PRIMER_AUTODETECTION_PRIMER_COL,
    "numreads": PRIMER_AUTODETECTION_NUMREADS_COL,
    "covbases": PRIMER_AUTODETECTION_COVBASES_COL,
    "coverage": PRIMER_AUTODETECTION_COVERAGE_COL,
}


def select_primer(input_path: Path, output_path: Path, sample_id: str) -> Tuple[pd.DataFrame, str]:
    """
    Concatenate the input coverage TSV files and select the record with the highest score.
    Return the selected primer data and the name of the selected primer.
    """
    # concatenate the coverage files to 1 single file with 1 single header
    # important: sort the primer records so that the latest version is at the top.
    # doing so, the latest version is used if primers have the same score
    primer_detection_df = concat(
        input_path=input_path,
        output_csv_path=output_path / f"{sample_id}{PRIMER_DETECTION_SUFFIX}",
        sortby_col=RNAME_COL,
        ascending=False,
        input_glob_pattern="*.tsv",
        input_sep="\t",
    )
    primer_detection_df.rename(columns=SAMTOOLS_COVERAGE_TO_PSGA_COLS, inplace=True)

    # extract the primer with the highest score
    # use deep copy to suppress the warning: SettingWithCopyWarning
    detected_primer_df_slice = primer_detection_df.loc[
        primer_detection_df[PRIMER_AUTODETECTION_PRIMER_SCORE_COL].idxmax()
    ].copy(deep=True)

    detected_primer = detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_COL]

    if detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_SCORE_COL] <= 0:
        # no primer can be used as no match was found.
        # overwrite the record
        detected_primer_df_slice[:] = 0
        detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_COL] = "none"
        # WARNING: this is a hack in order to run ncov.
        # certain samples could not have primers at all
        # (e.g. bam files in which adapters and primers were trimmed in the past)
        # ncov and artic expect a primer scheme and version as input parameters, though.
        # here we mock ncov by passing the latest ARTIC primer to make it run.
        # This won't have consequences because none of the sequences in this primer was found in the sample,
        # so they won't be trimmed
        # extract the first (=latest version) primer name matching ARTIC
        detected_primer = primer_detection_df[
            primer_detection_df[PRIMER_AUTODETECTION_PRIMER_COL].str.startswith("ARTIC_")
        ][PRIMER_AUTODETECTION_PRIMER_COL].iloc[0]

    return detected_primer_df_slice, detected_primer


def write_primer_data(
    output_path: Path,
    detected_primer_df_slice: pd.DataFrame,
    sample_id: str,
    primer_input: str,
) -> None:
    """
    Generate primer_data_df with the new columns and store it
    """
    primer_data_df = pd.DataFrame(detected_primer_df_slice).transpose().reset_index(drop=True)
    primer_data_df = primer_data_df.assign(
        **{
            PRIMER_AUTODETECTION_SAMPLE_ID_COL: sample_id,
            PRIMER_AUTODETECTION_PRIMER_INPUT_COL: primer_input,
        }
    )
    primer_data_df.to_csv(output_path / f"{sample_id}{PRIMER_DATA_SUFFIX}", index=False)


def write_selected_primer(output_path: Path, sample_id: str, detected_primer: str) -> None:
    """
    Store the primer scheme name/version
    """
    with open(output_path / f"{sample_id}_primer.txt", "w") as of:
        of.write(detected_primer)


def generate_primer_autodetection_output_files(
    input_path: Path, output_path: Path, sample_id: str, primer_input: str
) -> None:
    """
    Generate the primer autodetection output files
    """
    detected_primer_df_slice, detected_primer = select_primer(input_path, output_path, sample_id)
    write_primer_data(output_path, detected_primer_df_slice, sample_id, primer_input)
    write_selected_primer(output_path, sample_id, detected_primer)


@click.command()
@click.option(
    "--input-path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    required=True,
    help="input directory containing the csv files to concatenate",
)
@click.option(
    "--output-path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, writable=True),
    default=".",
    help="output directory to store the file to generate",
)
@click.option(
    "--sample-id",
    type=str,
    required=True,
    help="The sample id",
)
@click.option(
    "--primer-input",
    type=str,
    required=True,
    help="The primer name / version (e.g. ARTIC_V4) passed to the pipeline",
)
def primer_autodetection(
    input_path: str,
    output_path: str,
    sample_id: str,
    primer_input: str,
) -> None:
    """
    Generate the primer autodetection output files
    """
    generate_primer_autodetection_output_files(Path(input_path), Path(output_path), sample_id, primer_input)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    primer_autodetection()
