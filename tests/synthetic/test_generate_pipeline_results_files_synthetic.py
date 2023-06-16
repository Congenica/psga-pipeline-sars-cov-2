from pathlib import Path
import pytest
import json
from click.testing import CliRunner

from scripts.synthetic.generate_pipeline_results_files import generate_pipeline_results_files
from scripts.util.metadata import SAMPLE_ID
from tests.utils_tests import assert_csvs_are_equal, assert_jsons_are_equal


@pytest.mark.parametrize("pathogen", ["synthetic"])
@pytest.mark.parametrize(
    "use_s3_uri",
    [False, True],
)
@pytest.mark.parametrize(
    "metadata,exp_results_csv,exp_results_json,exp_resultfiles_json",
    [
        (
            "metadata_illumina.csv",
            "results_illumina.csv",
            "results_illumina.json",
            "resultfiles_illumina.json",
        ),
        (
            "metadata_ont.csv",
            "results_ont.csv",
            "results_ont.json",
            "resultfiles_ont.json",
        ),
        (
            "metadata_unknown.csv",
            "results_unknown.csv",
            "results_unknown.json",
            "resultfiles_unknown.json",
        ),
    ],
)
@pytest.mark.jira(identifier="1b81ce66-aa70-438e-986f-3684aaed652e", confirms="PSG-3621")
def test_generate_pipeline_results_files(
    tmp_path: Path,
    pipeline_results_files_data_path: Path,
    pathogen: str,
    metadata: str,
    exp_results_csv: str,
    exp_results_json: str,
    exp_resultfiles_json: str,
    use_s3_uri: bool,
):

    output_results_csv_file = tmp_path / "results.csv"
    output_results_json_file = tmp_path / "results.json"
    output_resultfiles_json_file = tmp_path / "resultfiles.json"

    output_path = "s3://bucket/analysis_run" if use_s3_uri else tmp_path

    args = [
        "--metadata-file",
        pipeline_results_files_data_path / pathogen / metadata,
        "--output-results-csv-file",
        output_results_csv_file,
        "--output-results-json-file",
        output_results_json_file,
        "--output-resultfiles-json-file",
        output_resultfiles_json_file,
        "--output-path",
        output_path,
    ]

    rv = CliRunner().invoke(
        generate_pipeline_results_files,
        args,
    )

    assert rv.exit_code == 0

    # assert results.csv
    exp_output_results_csv_file = pipeline_results_files_data_path / pathogen / exp_results_csv
    assert_csvs_are_equal(output_results_csv_file, exp_output_results_csv_file, SAMPLE_ID)

    # assert results.json
    exp_output_results_json_file = pipeline_results_files_data_path / pathogen / exp_results_json
    assert_jsons_are_equal(output_results_json_file, exp_output_results_json_file)

    # assert resultfiles.json
    with open(pipeline_results_files_data_path / pathogen / exp_resultfiles_json) as json_fd:
        exp_resultfiles_json_dict = json.load(json_fd)
    with open(output_resultfiles_json_file) as json_fd:
        calc_resultfiles_json_dict = json.load(json_fd)

    exp_resultfiles_json_dict_full_path = {
        sample_id: [{k: f"{output_path}/{v}" if k == "file" else v for k, v in d.items()} for d in list_of_dicts]
        for sample_id, list_of_dicts in exp_resultfiles_json_dict.items()
    }

    assert exp_resultfiles_json_dict_full_path == calc_resultfiles_json_dict
