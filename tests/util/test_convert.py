import pytest

from scripts.util.convert import csv_to_json
from tests.utils_tests import assert_jsons_are_equal


@pytest.mark.parametrize(
    "input_csv,expected_json,sample_id",
    [
        (
            "results_illumina_ont.csv",
            "results_illumina_ont.json",
            "SAMPLE_ID",
        ),
        (
            "results_unknown.csv",
            "results_unknown.json",
            "SAMPLE_ID",
        ),
    ],
)
def test_csv_to_json(tmp_path, test_data_path, input_csv, expected_json, sample_id):
    input_path = test_data_path / "pipeline_results_files" / input_csv
    expected_output_path = test_data_path / "pipeline_results_files" / expected_json
    output_path = tmp_path / "results.json"

    csv_to_json(input_path, output_path, sample_id)

    assert_jsons_are_equal(output_path, expected_output_path)
