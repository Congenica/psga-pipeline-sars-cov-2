from pathlib import Path
import pandas as pd
import click

from concat_csv import concat

SAMTOOLS_COVERAGE_TSV_EXTENSION = ".coverage.tsv"
PRIMER_DETECTION_SUFFIX = "_primer_detection.csv"
PRIMER_DATA_SUFFIX = "_primer_data.csv"
PRIMER_QC_PASS = "PASS"
PRIMER_QC_FAIL = "FAIL"

RNAME_COL = "#rname"
PRIMER_AUTODETECTION_SAMPLE_ID_COL = "sample_id"
PRIMER_AUTODETECTION_QC_COL = "primer_qc"
PRIMER_AUTODETECTION_PRIMER_INPUT_COL = "primer_input"
PRIMER_AUTODETECTION_PRIMER_COL = "primer_detected"
PRIMER_AUTODETECTION_NUMREADS_COL = "primer_numreads"
PRIMER_AUTODETECTION_COVBASES_COL = "primer_covbases"
PRIMER_AUTODETECTION_COVERAGE_COL = "primer_coverage"
EXPECTED_PRIMER_AUTODETECTION_HEADERS = {
    PRIMER_AUTODETECTION_SAMPLE_ID_COL,
    PRIMER_AUTODETECTION_QC_COL,
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
PRIMER_AUTODETECTION_PRIMER_SCORE_COL = PRIMER_AUTODETECTION_NUMREADS_COL  # primer_coverage might be more correct.
# columns to be renamed
SAMTOOLS_COVERAGE_TO_PSGA_COLS = {
    RNAME_COL: PRIMER_AUTODETECTION_PRIMER_COL,
    "numreads": PRIMER_AUTODETECTION_NUMREADS_COL,
    "covbases": PRIMER_AUTODETECTION_COVBASES_COL,
    "coverage": PRIMER_AUTODETECTION_COVERAGE_COL,
}


def select_primer(input_path: Path, output_path: Path, sample_id: str) -> pd.DataFrame:
    """
    Concatenate the input coverage TSV files and select the record with the highest score
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

    if detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_SCORE_COL] <= 0:
        # no primer can be used as no match was found.
        # overwrite the record
        detected_primer_df_slice[:] = 0
        detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_COL] = "none"

    return detected_primer_df_slice


def calculate_primer_qc(primer_input: str, primer_detected: str, primer_score: float) -> str:
    """
    Calculate the primer qc for the detected primer
    """
    if (primer_detected == "none" and primer_score > 0) or (primer_detected != "none" and primer_score == 0):
        raise ValueError(
            "Error in computing primer QC: "
            f"primer_detected={primer_detected}, primer_score={primer_score} are incompatible"
        )

    if primer_input == "unknown":
        primer_qc = PRIMER_QC_PASS if primer_score > 0 else PRIMER_QC_FAIL
    elif primer_input == "none":
        primer_qc = PRIMER_QC_FAIL if primer_score > 0 else PRIMER_QC_PASS
    else:
        primer_qc = PRIMER_QC_PASS if primer_input == primer_detected else PRIMER_QC_FAIL

    return primer_qc


def write_primer_data(
    output_path: Path,
    detected_primer_df_slice: pd.DataFrame,
    sample_id: str,
    primer_qc: str,
    primer_input: str,
) -> None:
    """
    Generate primer_data_df with the new columns and store it
    """
    primer_data_df = pd.DataFrame(detected_primer_df_slice).transpose().reset_index(drop=True)
    primer_data_df = primer_data_df.assign(
        **{
            PRIMER_AUTODETECTION_SAMPLE_ID_COL: sample_id,
            PRIMER_AUTODETECTION_QC_COL: primer_qc,
            PRIMER_AUTODETECTION_PRIMER_INPUT_COL: primer_input,
        }
    )
    primer_data_df.to_csv(output_path / f"{sample_id}{PRIMER_DATA_SUFFIX}", index=False)


def write_selected_primer(output_path: Path, sample_id: str, primer_input: str, primer: str, primer_qc: str) -> None:
    """
    Store the primer scheme name/version
    """
    with open(output_path / f"{sample_id}_primer_{primer_qc}.txt", "w") as of:
        if primer == "none" and primer_input == "none":
            # WARNING: this is a hack
            # certain samples could not have primers at all
            # (e.g. bam files in which adapters and primers were trimmed in the past)
            # ncov and artic expect a primer scheme and version as input parameters, though.
            # here we mock ncov by passing a primer scheme/version to make it run.
            # This won't have consequences because none of the sequences in this primer was found in the sample,
            # so they won't be trimmed
            of.write("ARTIC_V4")
        else:
            of.write(primer)


def generate_primer_autodetection_output_files(
    input_path: Path, output_path: Path, sample_id: str, primer_input: str
) -> None:
    """
    Generate the primer autodetection output files
    """
    detected_primer_df_slice = select_primer(input_path, output_path, sample_id)
    primer_score = detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_SCORE_COL]
    primer = detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_COL]

    primer_qc = calculate_primer_qc(primer_input, primer, primer_score)

    write_primer_data(output_path, detected_primer_df_slice, sample_id, primer_qc, primer_input)

    write_selected_primer(output_path, sample_id, primer_input, primer, primer_qc)


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
