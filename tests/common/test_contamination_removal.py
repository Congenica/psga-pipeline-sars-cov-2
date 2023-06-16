from pathlib import Path
import csv
import pytest
from click.testing import CliRunner

from scripts.common.contamination_removal import (
    get_contaminated_reads,
    write_rik_output_csv,
    process_rik,
    contamination_removal,
    CONTAMINATION_REMOVAL_SAMPLE_ID_COL,
    CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL,
)


def assert_rik_output_csv(output_path: Path, expected_rik_output_csv: dict[str, str]):
    with open(output_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        reader_list = list(reader)

    assert len(reader_list) == 1
    assert reader_list[0] == expected_rik_output_csv


@pytest.mark.parametrize(
    "input_file",
    [
        "rik_output_two_reads_exception.txt",
    ],
)
@pytest.mark.jira(identifier="16c0752b-146f-4ea1-ab63-4fc2c04d6575", confirms="PSG-3621")
def test_get_contaminated_reads_exception(contamination_removal_data_path: Path, input_file: str):
    input_path = contamination_removal_data_path / input_file
    with pytest.raises(ValueError, match="cannot be negative"):
        get_contaminated_reads(input_path)


@pytest.mark.parametrize(
    "input_file,expected_contaminated_reads",
    [
        ("rik_output_one_read.txt", 0),
        ("rik_output_two_reads.txt", 696),
    ],
)
@pytest.mark.jira(identifier="80d91164-7c98-4399-80ef-821627f8336f", confirms="PSG-3621")
def test_get_contaminated_reads(
    contamination_removal_data_path: Path, input_file: str, expected_contaminated_reads: int
):
    input_path = contamination_removal_data_path / input_file
    assert expected_contaminated_reads == get_contaminated_reads(input_path)


@pytest.mark.parametrize(
    "sample_id,contaminated_reads,expected_rik_output_csv",
    [
        (
            "a",
            10,
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: "a",
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: "10",
            },
        ),
        (
            "b",
            0,
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: "b",
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: "0",
            },
        ),
    ],
)
@pytest.mark.jira(identifier="12f810d0-8760-4e33-8912-db5dfff59675", confirms="PSG-3621")
def test_write_rik_output_csv(
    tmp_path: Path, sample_id: str, contaminated_reads: int, expected_rik_output_csv: dict[str, str]
):
    output_path = tmp_path / "rik_output.csv"
    write_rik_output_csv(output_path, sample_id, contaminated_reads)
    assert_rik_output_csv(output_path, expected_rik_output_csv)


@pytest.mark.parametrize(
    "sample_id,input_file,expected_rik_output_csv",
    [
        (
            "a",
            "rik_output_one_read.txt",
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: "a",
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: "0",
            },
        ),
        (
            "b",
            "rik_output_two_reads.txt",
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: "b",
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: "696",
            },
        ),
    ],
)
@pytest.mark.jira(identifier="15090941-15df-46d2-9a5d-5ffd660818d4", confirms="PSG-3621")
def test_process_rik(
    tmp_path: Path,
    contamination_removal_data_path: Path,
    sample_id: str,
    input_file: str,
    expected_rik_output_csv: dict[str, str],
):
    input_path = contamination_removal_data_path / input_file
    output_path = tmp_path / "rik_output.csv"
    process_rik(input_path, output_path, sample_id)
    assert_rik_output_csv(output_path, expected_rik_output_csv)


@pytest.mark.parametrize(
    "sample_id,input_file,expected_rik_output_csv",
    [
        (
            "a",
            "rik_output_one_read.txt",
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: "a",
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: "0",
            },
        ),
        (
            "b",
            "rik_output_two_reads.txt",
            {
                CONTAMINATION_REMOVAL_SAMPLE_ID_COL: "b",
                CONTAMINATION_REMOVAL_CONTAMINATED_READS_COL: "696",
            },
        ),
    ],
)
@pytest.mark.jira(identifier="5fb18309-ccdf-44ce-852d-d723d767b180", confirms="PSG-3621")
def test_contamination_removal(
    tmp_path: Path,
    contamination_removal_data_path: Path,
    sample_id: str,
    input_file: str,
    expected_rik_output_csv: dict[str, str],
):
    input_path = contamination_removal_data_path / input_file
    output_path = tmp_path / "rik_output.csv"
    rv = CliRunner().invoke(
        contamination_removal,
        [
            "--input-path",
            input_path,
            "--output-csv-path",
            output_path,
            "--sample-id",
            sample_id,
        ],
    )
    assert rv.exit_code == 0
    assert_rik_output_csv(output_path, expected_rik_output_csv)
