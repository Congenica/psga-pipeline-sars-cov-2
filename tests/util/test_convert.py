from pathlib import Path
import pytest

from scripts.util.convert import csv_to_json
from tests.utils_tests import assert_jsons_are_equal


@pytest.mark.parametrize("pathogen", ["sars_cov_2"])
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
@pytest.mark.jira(identifier="102bd1eb-18ec-4c15-b92e-f0ba55ccf61a", confirms="PSG-3621")
def test_csv_to_json(
    tmp_path: Path,
    pipeline_results_files_data_path: Path,
    pathogen: str,
    input_csv: str,
    expected_json: str,
    sample_id: str,
):
    input_path = pipeline_results_files_data_path / pathogen / input_csv
    expected_output_path = pipeline_results_files_data_path / pathogen / expected_json
    output_path = tmp_path / "results.json"

    csv_to_json(input_path, output_path, sample_id)

    assert_jsons_are_equal(output_path, expected_output_path)
