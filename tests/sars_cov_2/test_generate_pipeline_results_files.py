import pytest
import json
from click.testing import CliRunner

from scripts.sars_cov_2.generate_pipeline_results_files import generate_pipeline_results_files
from scripts.util.metadata import SAMPLE_ID
from tests.utils_tests import assert_csvs_are_equal, assert_jsons_are_equal


@pytest.mark.parametrize(
    "use_s3_uri",
    [False, True],
)
@pytest.mark.parametrize(
    "metadata,contamination_removal_csv,primer_autodetection_csv,ncov_qc_csv,ncov_typing_csv,pangolin_csv,"
    "analysis_run_name,sequencing_technology,exp_results_csv,exp_results_json,exp_resultfiles_json",
    [
        (
            "metadata_illumina.csv",
            "contamination_removal.csv",
            "primer_autodetection.csv",
            "ncov_test.qc.csv",
            "ncov_typing.csv",
            "all_lineages_report.csv",
            "just_a_name",
            "illumina",
            "results_illumina_ont.csv",
            "results_illumina_ont.json",
            "resultfiles_illumina.json",
        ),
        (
            "metadata_ont.csv",
            "contamination_removal.csv",
            "primer_autodetection.csv",
            "ncov_test.qc.csv",
            "ncov_typing.csv",
            "all_lineages_report.csv",
            "just_a_name",
            "ont",
            "results_illumina_ont.csv",
            "results_illumina_ont.json",
            "resultfiles_ont.json",
        ),
        (
            "metadata_unknown.csv",
            None,
            None,
            None,
            None,
            "all_lineages_report.csv",
            "just_a_name",
            "unknown",
            "results_unknown.csv",
            "results_unknown.json",
            "resultfiles_unknown.json",
        ),
    ],
)
def test_generate_pipeline_results_files(
    tmp_path,
    test_data_path,
    metadata,
    contamination_removal_csv,
    primer_autodetection_csv,
    ncov_qc_csv,
    ncov_typing_csv,
    pangolin_csv,
    analysis_run_name,
    sequencing_technology,
    exp_results_csv,
    exp_results_json,
    exp_resultfiles_json,
    use_s3_uri,
):

    output_results_csv_file = tmp_path / "results.csv"
    output_results_json_file = tmp_path / "results.json"
    output_resultfiles_json_file = tmp_path / "resultfiles.json"

    output_path = "s3://bucket/analysis_run" if use_s3_uri else tmp_path

    args = [
        "--analysis-run-name",
        analysis_run_name,
        "--metadata-file",
        test_data_path / "pipeline_results_files" / metadata,
        "--pangolin-csv-file",
        test_data_path / "pipeline_results_files" / pangolin_csv,
        "--output-results-csv-file",
        output_results_csv_file,
        "--output-results-json-file",
        output_results_json_file,
        "--output-resultfiles-json-file",
        output_resultfiles_json_file,
        "--output-path",
        output_path,
        "--sequencing-technology",
        sequencing_technology,
    ]

    if ncov_qc_csv:
        args.extend(
            [
                "--contamination-removal-csv-file",
                test_data_path / "pipeline_results_files" / contamination_removal_csv,
                "--primer-autodetection-csv-file",
                test_data_path / "pipeline_results_files" / primer_autodetection_csv,
                "--ncov-qc-csv-file",
                test_data_path / "pipeline_results_files" / ncov_qc_csv,
                "--ncov-typing-csv-file",
                test_data_path / "pipeline_results_files" / ncov_typing_csv,
            ]
        )
    rv = CliRunner().invoke(
        generate_pipeline_results_files,
        args,
    )

    assert rv.exit_code == 0

    # assert results.csv
    exp_output_results_csv_file = test_data_path / "pipeline_results_files" / exp_results_csv
    assert_csvs_are_equal(output_results_csv_file, exp_output_results_csv_file, SAMPLE_ID)

    # assert results.json
    exp_output_results_json_file = test_data_path / "pipeline_results_files" / exp_results_json
    assert_jsons_are_equal(output_results_json_file, exp_output_results_json_file)

    # assert resultfiles.json
    with open(test_data_path / "pipeline_results_files" / exp_resultfiles_json) as json_fd:
        exp_resultfiles_json_dict = json.load(json_fd)
    with open(output_resultfiles_json_file) as json_fd:
        calc_resultfiles_json_dict = json.load(json_fd)

    exp_resultfiles_json_dict_full_path = {
        sample_id: [{k: f"{output_path}/{v}" if k == "file" else v for k, v in d.items()} for d in list_of_dicts]
        for sample_id, list_of_dicts in exp_resultfiles_json_dict.items()
    }

    assert exp_resultfiles_json_dict_full_path == calc_resultfiles_json_dict
