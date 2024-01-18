from pathlib import Path
from typing import Any
import csv
import pickle
import gzip
from dataclasses import dataclass, field
import pandas as pd
import click
from Bio import SeqIO
import ahocorasick

from concat_csv import concat

from scripts.primer_colsprimer_cols import (
    PRIMER_INDEX_COLS,
    PRIMER_NAME,
    PICKLE_PATH,
    TOTAL_NUM_PRIMER,
    PRIMER_AUTODETECTION_SAMPLE_ID_COL,
    PRIMER_AUTODETECTION_PRIMER_INPUT_COL,
    PRIMER_AUTODETECTION_PRIMER_COL,
    PRIMER_AUTODETECTION_NUMREADS_COL,
    PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL,
    PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL,
    PRIMER_AUTODETECTION_COVERAGE_COL,
)

UNKNOWN = "unknown"

COVERAGE_SUFFIX = ".coverage.csv"
PRIMER_DETECTION_SUFFIX = "_primer_detection.csv"
PRIMER_DATA_SUFFIX = "_primer_data.csv"

PRIMER_AUTODETECTION_PRIMER_SCORE_COL = PRIMER_AUTODETECTION_NUMREADS_COL


@dataclass
class PrimerAutomaton:
    data: dict = field(metadata={"required": True})
    automaton: ahocorasick.Automaton = field(metadata={"required": True})


def load_pickle(input_path: Path) -> Any:
    """
    Load a Python object using pickle
    """
    with open(input_path, "rb") as pickle_in:
        return pickle.load(pickle_in)


def build_primers_automaton(primer_index: Path) -> dict[str, PrimerAutomaton]:
    """
    Build the primers automaton dictionary
    """
    primers_automaton = {}

    with open(primer_index, "r", newline="") as index_csv:
        reader = csv.DictReader(index_csv, fieldnames=PRIMER_INDEX_COLS)
        # the first row is the header
        _ = next(reader, None)
        for row in reader:
            primer = row[PRIMER_NAME]
            primer_pickle_path = Path(row[PICKLE_PATH])
            data = {
                PRIMER_AUTODETECTION_PRIMER_COL: primer,
                TOTAL_NUM_PRIMER: int(row[TOTAL_NUM_PRIMER]),
                PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 0,
            }

            # this can be optimised by generating the automaton at build time (one-off)
            automaton = load_pickle(primer_pickle_path)
            primers_automaton[primer] = PrimerAutomaton(
                data=data,
                automaton=automaton,
            )

    return primers_automaton


def count_primer_matches(sample_fastq: Path, primers_automaton: dict) -> dict:
    """
    Perform an exact search of the primer sequences (all schemes) in the sample_fastq
    """
    # Structure for storing the unique primer sequences
    unique_hits: dict[str, set] = {p: set() for p in primers_automaton}

    # look up the primer sequences in the automaton
    # reading the sample fastq once for all
    with gzip.open(sample_fastq, "rt") as fastq_file:
        # for each sample read, it counts the number of primers found for each primer scheme.
        # This search is performed using an ahocorasick structure and optimised to search
        # at the beginning of the read only
        for sample_sequence in SeqIO.parse(fastq_file, "fastq"):
            sample_read = str(sample_sequence.seq)

            for primer in primers_automaton:
                primer_auto = primers_automaton[primer]
                # search for primers matching the prefix of sample_read
                # (the wildcard argument is mandatory, but unused - no base will ever be '?'), e.g.
                # keys: ['AGA', 'AAGA', 'ATGA', 'GTAT']
                # sample_read: 'AAGATT'
                # output: ['AAGA']
                found_primer = set(primer_auto.automaton.keys(sample_read, "?", ahocorasick.MATCH_AT_MOST_PREFIX))
                # search for primers matching the prefix of sample_read (wildcard 'N' is used).
                # This enables the counting of primers upon uncertainty in sample reads. e.g.
                # keys: ['AGA', 'AAGA', 'ATGA', 'GTAT']
                # sample_read: 'ANGATT'
                # output: ['AAGA', 'ATGA']
                # Depending on the number of Ns in sample_read, this can return:
                # - "sample_read does not contain N in primer sequence"
                #     => found_primer = {} or {primer}, ambiguous_primers = {}
                # - "sample_read contains Ns in primer sequence"
                #     => found_primer = {}, ambiguous_primers = {*}
                # The found primer (if exists) is removed as that's an exact match
                ambiguous_primers = (
                    set(primer_auto.automaton.keys(sample_read, "N", ahocorasick.MATCH_AT_MOST_PREFIX)) - found_primer
                )

                # increment the read count as necessary
                if found_primer:
                    primer_auto.data[PRIMER_AUTODETECTION_NUMREADS_COL] += 1
                    unique_hits[primer].update(found_primer)

                if ambiguous_primers:
                    primer_auto.data[PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL] += 1

    # Add the unique number of primers found in the sample reads
    for primer in unique_hits:
        primers_automaton[primer].data[PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL] = len(unique_hits[primer])

    return unique_hits


def compute_primer_data(sample_fastq: Path, primers_automaton: dict) -> dict:
    """
    Perform an exact search of the primer sequences (all schemes) in the sample_fastq
    """
    count_primer_matches(sample_fastq, primers_automaton)

    return primers_automaton


def generate_metrics(primer_index: Path, sample_fastq: Path, output_path: Path) -> None:
    """
    Generate metrics for all primers.
    Metrics are saved to separate CSV files
    """
    primer_coverage_csv_fieldnames = [
        PRIMER_AUTODETECTION_PRIMER_COL,
        TOTAL_NUM_PRIMER,
        PRIMER_AUTODETECTION_NUMREADS_COL,
        PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL,
        PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL,
    ]

    primers_automaton = build_primers_automaton(primer_index)

    primers_automaton = compute_primer_data(sample_fastq, primers_automaton)

    for primer in primers_automaton:
        with open(output_path / f"{primer}{COVERAGE_SUFFIX}", "w", newline="") as primer_coverage_csv:
            writer = csv.DictWriter(primer_coverage_csv, fieldnames=primer_coverage_csv_fieldnames)
            writer.writeheader()
            writer.writerow(primers_automaton[primer].data)


def select_primer(output_path: Path, sample_id: str, primer_input: str) -> tuple[pd.DataFrame, str]:
    """
    Concatenate the input coverage files and select the record with the highest score.
    Return the selected primer data and the name of the selected primer.
    """
    # concatenate the coverage files to 1 single file with 1 single header
    # important: sort the primer records so that the latest version is at the top.
    # doing so, the latest version is used if primers have the same score
    primer_detection_df = concat(
        input_path=output_path,
        output_csv_path=output_path / f"{sample_id}{PRIMER_DETECTION_SUFFIX}",
        sortby_col=PRIMER_AUTODETECTION_PRIMER_COL,
        ascending=False,
        input_glob_pattern=f"*{COVERAGE_SUFFIX}",
        input_sep=",",
    )

    # calculate coverage as a percentage of primer sequences hit for a certain primer scheme
    primer_detection_df[PRIMER_AUTODETECTION_COVERAGE_COL] = (
        primer_detection_df[PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL] / primer_detection_df[TOTAL_NUM_PRIMER]
    )
    primer_detection_df.drop(columns=[TOTAL_NUM_PRIMER], inplace=True)

    # extract the primer with the highest score
    # use deep copy to suppress the warning: settingWithCopyWarning
    detected_primer_df_slice = primer_detection_df.loc[
        primer_detection_df[PRIMER_AUTODETECTION_PRIMER_SCORE_COL].idxmax()
    ].copy(deep=True)

    detected_primer = detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_COL]

    if detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_SCORE_COL] <= 0:
        # no primer can be used as no match was found.
        # overwrite the record
        detected_primer_df_slice[:] = 0
        detected_primer_df_slice[PRIMER_AUTODETECTION_PRIMER_COL] = "none"
        # WARNING: this is a hack in order to run ncov, in the case that the input primer is unknown.
        # certain samples could not have primers at all
        # (e.g. bam files in which adapters and primers were trimmed in the past)
        # ncov and artic expect a primer scheme and version as input parameters, though.
        # here we mock ncov passing the latest ARTIC primer to make it run.
        # This won't have consequences because no primer sequences were found in the sample,
        # so they won't be trimmed.
        # extract the first (=latest version) primer name matching ARTIC
        detected_primer = primer_detection_df[
            primer_detection_df[PRIMER_AUTODETECTION_PRIMER_COL].str.startswith("ARTIC_")
        ][PRIMER_AUTODETECTION_PRIMER_COL].iloc[0]

    # if known, the input primer has precedence over this primer autodetection.
    # both the input and the detected primers are reported, so that the infrastructure can raise a QC warning
    selected_primer = detected_primer if primer_input == UNKNOWN else primer_input

    return detected_primer_df_slice, selected_primer


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


def write_selected_primer(output_path: Path, sample_id: str, selected_primer: str) -> None:
    """
    Store the primer scheme name/version
    """
    with open(output_path / f"{sample_id}_primer.txt", "w") as of:
        of.write(selected_primer)


def generate_primer_autodetection_output_files(output_path: Path, sample_id: str, primer_input: str) -> None:
    """
    Generate the primer autodetection output files
    """
    detected_primer_df_slice, selected_primer = select_primer(output_path, sample_id, primer_input)
    write_primer_data(output_path, detected_primer_df_slice, sample_id, primer_input)
    write_selected_primer(output_path, sample_id, selected_primer)


@click.command()
@click.option(
    "--primer-index",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="path to the primer index file",
)
@click.option(
    "--sample-fastq",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="path to the primer index file",
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
    primer_index: str,
    sample_fastq: str,
    output_path: str,
    sample_id: str,
    primer_input: str,
) -> None:
    """
    Generate the primer autodetection output files
    """
    output_path_obj = Path(output_path)
    generate_metrics(Path(primer_index), Path(sample_fastq), output_path_obj)
    generate_primer_autodetection_output_files(output_path_obj, sample_id, primer_input)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    primer_autodetection()
