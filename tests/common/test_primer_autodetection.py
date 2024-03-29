from pathlib import Path
import shutil
import csv
import pytest

from click.testing import CliRunner

from app.scripts.fetch_primers import (
    create_automaton,
)

from app.scripts.primer_autodetection import (
    load_pickle,
    build_primers_automaton,
    count_primer_matches,
    compute_primer_data,
    generate_metrics,
    select_primer,
    write_primer_data,
    write_selected_primer,
    generate_primer_autodetection_output_files,
    primer_autodetection,
    PrimerAutomaton,
    PRIMER_AUTODETECTION_PRIMER_SCORE_COL,
    COVERAGE_SUFFIX,
    PRIMER_DETECTION_SUFFIX,
    PRIMER_DATA_SUFFIX,
    UNKNOWN,
)
from app.scripts.fetch_primers import (
    SARS_COV_2,
    SCHEME,
    PICKLE,
)
from app.scripts.primer_cols import (
    PRIMER_INDEX_COLS,
    PRIMER_NAME,
    FASTA_PATH,
    PICKLE_PATH,
    TOTAL_NUM_PRIMER,
    PRIMER_AUTODETECTION_SAMPLE_ID_COL,
    PRIMER_AUTODETECTION_PRIMER_COL,
    PRIMER_AUTODETECTION_NUMREADS_COL,
    PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL,
    PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL,
)
from tests.utils_tests import assert_csvs_are_equal, assert_files_are_equal

DEFAULT_PRIMER = "ARTIC_V4-1"
PRIMER_TEST = "primer_test"


def copy_with_wildcard(input_path: Path, output_path: Path, suffix: str):
    for my_path in input_path.glob(f"*{suffix}"):
        shutil.copy(my_path, output_path)


def prefix_path_to_index(input_path: Path, output_path: Path, path_prefix: Path):
    primers = []
    with open(input_path, "r", newline="") as fin, open(output_path, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=PRIMER_INDEX_COLS)
        writer.writeheader()
        for row in reader:
            row[FASTA_PATH] = f"{path_prefix}{row[FASTA_PATH]}"
            row[PICKLE_PATH] = f"{path_prefix}{row[PICKLE_PATH]}"
            writer.writerow(row)
            primers.append(row[PRIMER_NAME])
    return primers


def assert_primer_coverage(primer: str, tmp_path: Path, input_path: Path):
    output_file = f"{primer}{COVERAGE_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_csvs_are_equal(output_path, expected_output_path, PRIMER_AUTODETECTION_PRIMER_COL)


def assert_primer_detection(sample_id: str, tmp_path: Path, input_path: Path):
    output_file = f"{sample_id}{PRIMER_DETECTION_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_csvs_are_equal(output_path, expected_output_path, PRIMER_AUTODETECTION_PRIMER_COL)


def assert_primer_data(sample_id: str, tmp_path: Path, input_path: Path):
    output_file = f"{sample_id}{PRIMER_DATA_SUFFIX}"
    output_path = tmp_path / output_file
    expected_output_path = input_path / output_file
    assert_csvs_are_equal(output_path, expected_output_path, PRIMER_AUTODETECTION_SAMPLE_ID_COL)


def assert_selected_primer_file(sample_id: str, primer_autodetection_data_path: Path, tmp_path: Path, found_dir: str):
    primer_file = f"{sample_id}_primer.txt"
    expected_output_path = primer_autodetection_data_path / found_dir / primer_file
    output_path = tmp_path / primer_file
    assert_files_are_equal(output_path, expected_output_path)


@pytest.fixture
def primers_automaton_fixture(primer_autodetection_primer_schemes_data_path: Path) -> dict[str, PrimerAutomaton]:
    automaton_dictionary = {
        "ARTIC_V4": PrimerAutomaton(
            data={
                PRIMER_AUTODETECTION_PRIMER_COL: "ARTIC_V4",
                TOTAL_NUM_PRIMER: 4,
                PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 0,
            },
            automaton=load_pickle(
                primer_autodetection_primer_schemes_data_path
                / "ARTIC"
                / SARS_COV_2
                / "V4"
                / f"{SARS_COV_2}.{SCHEME}.{PICKLE}"
            ),
        ),
        "Midnight-ONT_V2": PrimerAutomaton(
            data={
                PRIMER_AUTODETECTION_PRIMER_COL: "Midnight-ONT_V2",
                TOTAL_NUM_PRIMER: 4,
                PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 0,
            },
            automaton=load_pickle(
                primer_autodetection_primer_schemes_data_path
                / "Midnight-ONT"
                / SARS_COV_2
                / "V2"
                / f"{SARS_COV_2}.{SCHEME}.{PICKLE}"
            ),
        ),
    }
    return automaton_dictionary


def assert_primer_automaton(
    primers_automaton: dict[str, PrimerAutomaton], expected_primers_automaton: dict[str, PrimerAutomaton]
):
    primers = list(primers_automaton.keys())
    expected_primers = list(expected_primers_automaton.keys())
    assert primers == expected_primers

    def _get_automaton_items(automaton_obj):
        return [t for t in automaton_obj.items()]

    for primer in primers:
        # compare the automaton items (primer sequences)
        automaton_items = _get_automaton_items(primers_automaton[primer].automaton)
        expected_automaton_items = _get_automaton_items(expected_primers_automaton[primer].automaton)
        assert automaton_items == expected_automaton_items
        # compare the rest of the data
        assert primers_automaton[primer].data == expected_primers_automaton[primer].data


@pytest.mark.parametrize(
    "index",
    [
        f"{SARS_COV_2}_primer_index.csv",
    ],
)
@pytest.mark.jira(identifier="7ea6a401-37f4-4e7d-8c99-b1b4276adfa1", confirms="PSG-3621")
def test_build_primers_automaton(
    tmp_path: Path,
    primer_autodetection_primer_schemes_data_path: Path,
    primer_autodetection_data_path: Path,
    index: str,
    primers_automaton_fixture: dict[str, PrimerAutomaton],
):
    orig_index_path = primer_autodetection_primer_schemes_data_path / index
    tmp_index_path = tmp_path / index
    prefix_path_to_index(orig_index_path, tmp_index_path, primer_autodetection_data_path)

    primers_automaton = build_primers_automaton(tmp_index_path)

    assert_primer_automaton(primers_automaton, primers_automaton_fixture)


@pytest.mark.parametrize(
    "sample,fasta,expected_unique_hits,expected_data",
    [
        (
            "sample_counting.fastq.gz",
            "sample_counting_primers.fasta",
            {PRIMER_TEST: {"AAGA"}},
            {
                PRIMER_AUTODETECTION_NUMREADS_COL: 2,
                PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 1,
                PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 1,
            },
        )
    ],
)
@pytest.mark.jira(identifier="5d595f04-ef60-4c98-ae45-14000dbb4143", confirms="PSG-3621")
def test_count_primer_matches(
    primer_autodetection_sample_dir_data_path: Path,
    sample: str,
    fasta: str,
    expected_unique_hits: dict[str, set],
    expected_data: dict[str, int],
):
    # build a simplified primer_automaton obj for this test
    primers_automaton = {
        PRIMER_TEST: PrimerAutomaton(
            data={
                PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 0,
            },
            automaton=create_automaton(primer_autodetection_sample_dir_data_path / fasta),
        ),
    }
    unique_hits = count_primer_matches(primer_autodetection_sample_dir_data_path / sample, primers_automaton)

    assert unique_hits == expected_unique_hits
    assert primers_automaton[PRIMER_TEST].data == expected_data


@pytest.mark.parametrize(
    "sample_id,expected_data",
    [
        (
            # ARTIC_V4 SAMPLE
            "9729bce7-f0a9-4617-b6e0-6145307741d1",
            {
                "ARTIC_V4": {
                    PRIMER_AUTODETECTION_PRIMER_COL: "ARTIC_V4",
                    TOTAL_NUM_PRIMER: 4,
                    PRIMER_AUTODETECTION_NUMREADS_COL: 3,
                    PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 2,
                    PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 1,
                },
                "Midnight-ONT_V2": {
                    PRIMER_AUTODETECTION_PRIMER_COL: "Midnight-ONT_V2",
                    TOTAL_NUM_PRIMER: 4,
                    PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                    PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                    PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 0,
                },
            },
        ),
        (
            # SAMPLE without actual primers
            "a0446f6f-7d24-478c-8d92-7c77036930d8",
            {
                "ARTIC_V4": {
                    PRIMER_AUTODETECTION_PRIMER_COL: "ARTIC_V4",
                    TOTAL_NUM_PRIMER: 4,
                    PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                    PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                    PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 1,
                },
                "Midnight-ONT_V2": {
                    PRIMER_AUTODETECTION_PRIMER_COL: "Midnight-ONT_V2",
                    TOTAL_NUM_PRIMER: 4,
                    PRIMER_AUTODETECTION_NUMREADS_COL: 0,
                    PRIMER_AUTODETECTION_UNIQUE_NUMREADS_COL: 0,
                    PRIMER_AUTODETECTION_AMBIGUOUS_NUMREADS_COL: 0,
                },
            },
        ),
    ],
)
@pytest.mark.jira(identifier="2116ec3e-87c2-40e3-bd70-aa369d4d79e2", confirms="PSG-3621")
def test_compute_primer_data(
    primer_autodetection_sample_dir_data_path: Path,
    primers_automaton_fixture: dict[str, PrimerAutomaton],
    sample_id: str,
    expected_data: dict[str, dict],
):
    sample_fastq = primer_autodetection_sample_dir_data_path / f"{sample_id}.fastq.gz"
    primers_automaton = compute_primer_data(sample_fastq, primers_automaton_fixture)
    for primer in primers_automaton:
        assert primers_automaton[primer].data == expected_data[primer]


@pytest.mark.parametrize(
    "index,sample_id",
    [
        (
            f"{SARS_COV_2}_primer_index.csv",
            # ARTIC_V4 SAMPLE
            "9729bce7-f0a9-4617-b6e0-6145307741d1",
        ),
        (
            f"{SARS_COV_2}_primer_index.csv",
            # SAMPLE without actual primers
            "a0446f6f-7d24-478c-8d92-7c77036930d8",
        ),
    ],
)
@pytest.mark.jira(identifier="61bf55d8-2717-44ce-bea9-f9ef6c05ee73", confirms="PSG-3621")
def test_generate_metrics(
    tmp_path: Path,
    primer_autodetection_data_path: Path,
    primer_autodetection_sample_dir_data_path: Path,
    primer_autodetection_primer_schemes_data_path: Path,
    index: str,
    sample_id: str,
):
    expected_output_path = primer_autodetection_sample_dir_data_path / "expected_output" / sample_id
    orig_index_path = primer_autodetection_primer_schemes_data_path / index
    tmp_index_path = tmp_path / index

    primers = prefix_path_to_index(orig_index_path, tmp_index_path, primer_autodetection_data_path)
    sample_fastq = primer_autodetection_sample_dir_data_path / f"{sample_id}.fastq.gz"

    generate_metrics(tmp_index_path, sample_fastq, tmp_path)

    for primer in primers:
        assert_primer_coverage(primer, tmp_path, expected_output_path)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input,expected_score,expected_detected_primer,expected_selected_primer",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "ARTIC_V3", 423, "ARTIC_V4-1", "ARTIC_V3"),
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", UNKNOWN, 423, "ARTIC_V4-1", "ARTIC_V4-1"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "ARTIC_V3", 0, "none", "ARTIC_V3"),
        # test the hack
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", UNKNOWN, 0, "none", DEFAULT_PRIMER),
    ],
)
@pytest.mark.jira(identifier="f370404a-c340-4677-bc83-507657db0d80", confirms="PSG-3621")
def test_select_primer(
    tmp_path: Path,
    primer_autodetection_data_path: Path,
    found_dir: str,
    sample_id: str,
    primer_input: str,
    expected_score: int,
    expected_detected_primer: str,
    expected_selected_primer: str,
):
    input_path = primer_autodetection_data_path / found_dir
    copy_with_wildcard(input_path, tmp_path, COVERAGE_SUFFIX)
    detected_primer_df, selected_primer = select_primer(tmp_path, sample_id, primer_input)

    assert expected_detected_primer == detected_primer_df[PRIMER_AUTODETECTION_PRIMER_COL]
    assert expected_score == detected_primer_df[PRIMER_AUTODETECTION_PRIMER_SCORE_COL]
    assert expected_selected_primer == selected_primer
    assert_primer_detection(sample_id, tmp_path, input_path)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown"),
    ],
)
@pytest.mark.jira(identifier="a52a1a34-0d96-405c-8782-4cdfcf7be06e", confirms="PSG-3621")
def test_write_primer_data(
    tmp_path: Path,
    primer_autodetection_data_path: Path,
    found_dir: str,
    sample_id: str,
    primer_input: str,
):
    input_path = primer_autodetection_data_path / found_dir
    copy_with_wildcard(input_path, tmp_path, COVERAGE_SUFFIX)
    detected_primer_df, _ = select_primer(tmp_path, sample_id, primer_input)
    write_primer_data(tmp_path, detected_primer_df, sample_id, primer_input)
    assert_primer_data(sample_id, tmp_path, input_path)


@pytest.mark.parametrize(
    "found_dir,sample_id,detected_primer",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "ARTIC_V4-1"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", DEFAULT_PRIMER),
    ],
)
@pytest.mark.jira(identifier="6831b9e0-7803-419f-8e08-0d554d1f6cf5", confirms="PSG-3621")
def test_write_selected_primer(
    tmp_path: Path,
    primer_autodetection_data_path: Path,
    found_dir: str,
    sample_id: str,
    detected_primer: str,
):
    write_selected_primer(tmp_path, sample_id, detected_primer)
    assert_selected_primer_file(sample_id, primer_autodetection_data_path, tmp_path, found_dir)


@pytest.mark.parametrize(
    "found_dir,sample_id,primer_input",
    [
        ("found", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown"),
        ("not_found", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown"),
    ],
)
@pytest.mark.jira(identifier="7d937435-6fda-49dc-80c2-1e74142a8d65", confirms="PSG-3621")
def test_generate_primer_autodetection_output_files(
    tmp_path: Path,
    primer_autodetection_data_path: Path,
    found_dir: str,
    sample_id: str,
    primer_input: str,
):
    input_path = primer_autodetection_data_path / found_dir
    copy_with_wildcard(input_path, tmp_path, COVERAGE_SUFFIX)
    generate_primer_autodetection_output_files(tmp_path, sample_id, primer_input)
    assert_primer_detection(sample_id, tmp_path, input_path)
    assert_primer_data(sample_id, tmp_path, input_path)
    assert_selected_primer_file(sample_id, primer_autodetection_data_path, tmp_path, found_dir)


@pytest.mark.parametrize(
    "index,sample_id,primer_input",
    [
        (f"{SARS_COV_2}_primer_index.csv", "9729bce7-f0a9-4617-b6e0-6145307741d1", "unknown"),
        (f"{SARS_COV_2}_primer_index.csv", "a0446f6f-7d24-478c-8d92-7c77036930d8", "unknown"),
    ],
)
@pytest.mark.jira(identifier="e4021c9e-609f-45b5-b366-9943b1146148", confirms="PSG-3621")
def test_primer_autodetection(
    tmp_path: Path,
    primer_autodetection_data_path: Path,
    primer_autodetection_sample_dir_data_path: Path,
    primer_autodetection_primer_schemes_data_path: Path,
    samples_dir: str,
    index: str,
    sample_id: str,
    primer_input: str,
):
    expected_output_path = primer_autodetection_sample_dir_data_path / "expected_output" / sample_id
    orig_index_path = primer_autodetection_primer_schemes_data_path / index
    tmp_index_path = tmp_path / index
    sample_fastq = primer_autodetection_sample_dir_data_path / f"{sample_id}.fastq.gz"

    prefix_path_to_index(orig_index_path, tmp_index_path, primer_autodetection_data_path)

    rv = CliRunner().invoke(
        primer_autodetection,
        [
            "--primer-index",
            tmp_index_path,
            "--sample-fastq",
            sample_fastq,
            "--output-path",
            tmp_path,
            "--sample-id",
            sample_id,
            "--primer-input",
            primer_input,
        ],
    )
    assert rv.exit_code == 0

    assert_primer_detection(sample_id, tmp_path, expected_output_path)
    assert_primer_data(sample_id, tmp_path, expected_output_path)
    assert_selected_primer_file(
        sample_id, primer_autodetection_data_path, tmp_path, f"{samples_dir}/expected_output/{sample_id}"
    )
