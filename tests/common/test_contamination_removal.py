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


def assert_rik_output_csv(output_path, expected_rik_output_csv):
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
def test_get_contaminated_reads_exception(contamination_removal_data_path, input_file):
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
def test_get_contaminated_reads(contamination_removal_data_path, input_file, expected_contaminated_reads):
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
def test_write_rik_output_csv(tmp_path, sample_id, contaminated_reads, expected_rik_output_csv):
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
def test_process_rik(tmp_path, contamination_removal_data_path, sample_id, input_file, expected_rik_output_csv):
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
def test_contamination_removal(
    tmp_path, contamination_removal_data_path, sample_id, input_file, expected_rik_output_csv
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
