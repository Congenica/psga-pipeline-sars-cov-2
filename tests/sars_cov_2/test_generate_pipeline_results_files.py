from pathlib import Path
import pytest
import json
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.sars_cov_2.generate_pipeline_results_files import (
    generate_pipeline_results_files,
    SAMPLES_UNKNOWN_NCOV_QC_FILE,
    SAMPLES_FAILED_NCOV_QC_FILE,
    SAMPLES_PASSED_NCOV_QC_FILE,
    SAMPLES_UNKNOWN_PANGOLIN_FILE,
    SAMPLES_FAILED_PANGOLIN_FILE,
    SAMPLES_PASSED_PANGOLIN_FILE,
    UNKNOWN_NCOV,
    FAILED_NCOV,
    PASSED_NCOV,
    UNKNOWN_PANGOLIN,
    FAILED_PANGOLIN,
    PASSED_PANGOLIN,
)
from scripts.util.metadata import SAMPLE_ID
from utils_tests import read_samples_from_file


def check_sample_list(input_path, expected_samples):
    processed_samples = read_samples_from_file(input_path)
    assert sorted(expected_samples) == sorted(processed_samples)


@pytest.mark.parametrize(
    "metadata,ncov_csv,pangolin_csv,analysis_run_name,sequencing_technology,"
    "exp_lists,exp_results_csv,exp_resultfiles_json",
    [
        (
            "metadata_illumina.csv",
            "ncov_test.qc.csv",
            "all_lineages_report.csv",
            "just_a_name",
            "illumina",
            {
                UNKNOWN_NCOV: ["985347c5-ff6a-454c-ac34-bc353d05dd70"],
                FAILED_NCOV: ["0774181d-fb20-4a73-b887-38af7eda9b38"],
                PASSED_NCOV: [
                    "e80f3c63-d139-4c14-bd72-7a43897ab40d",
                    "56a63f60-764d-4fc7-8764-a46023cbe324",
                    "a0951432-cd94-45b5-96d5-b721c037a451",
                ],
                UNKNOWN_PANGOLIN: ["a0951432-cd94-45b5-96d5-b721c037a451"],
                FAILED_PANGOLIN: ["e80f3c63-d139-4c14-bd72-7a43897ab40d"],
                PASSED_PANGOLIN: ["56a63f60-764d-4fc7-8764-a46023cbe324"],
            },
            "results_illumina_ont.csv",
            "resultfiles_illumina.json",
        ),
        (
            "metadata_ont.csv",
            "ncov_test.qc.csv",
            "all_lineages_report.csv",
            "just_a_name",
            "ont",
            {
                UNKNOWN_NCOV: ["985347c5-ff6a-454c-ac34-bc353d05dd70"],
                FAILED_NCOV: ["0774181d-fb20-4a73-b887-38af7eda9b38"],
                PASSED_NCOV: [
                    "e80f3c63-d139-4c14-bd72-7a43897ab40d",
                    "56a63f60-764d-4fc7-8764-a46023cbe324",
                    "a0951432-cd94-45b5-96d5-b721c037a451",
                ],
                UNKNOWN_PANGOLIN: ["a0951432-cd94-45b5-96d5-b721c037a451"],
                FAILED_PANGOLIN: ["e80f3c63-d139-4c14-bd72-7a43897ab40d"],
                PASSED_PANGOLIN: ["56a63f60-764d-4fc7-8764-a46023cbe324"],
            },
            "results_illumina_ont.csv",
            "resultfiles_ont.json",
        ),
        (
            "metadata_unknown.csv",
            None,
            "all_lineages_report.csv",
            "just_a_name",
            "unknown",
            {
                UNKNOWN_PANGOLIN: [
                    "985347c5-ff6a-454c-ac34-bc353d05dd70",
                    "0774181d-fb20-4a73-b887-38af7eda9b38",
                    "a0951432-cd94-45b5-96d5-b721c037a451",
                ],
                FAILED_PANGOLIN: ["e80f3c63-d139-4c14-bd72-7a43897ab40d"],
                PASSED_PANGOLIN: ["56a63f60-764d-4fc7-8764-a46023cbe324"],
            },
            "results_unknown.csv",
            "resultfiles_unknown.json",
        ),
    ],
)
def test_generate_pipeline_results_files(
    tmp_path,
    test_data_path,
    metadata,
    ncov_csv,
    pangolin_csv,
    analysis_run_name,
    sequencing_technology,
    exp_lists,
    exp_results_csv,
    exp_resultfiles_json,
):

    output_csv_file = Path(tmp_path / "results.csv")
    output_json_file = Path(tmp_path / "resultfiles.json")

    args = [
        "--analysis-run-name",
        analysis_run_name,
        "--metadata-file",
        test_data_path / "pipeline_results_files" / metadata,
        "--pangolin-csv-file",
        test_data_path / "pipeline_results_files" / pangolin_csv,
        "--output-csv-file",
        output_csv_file,
        "--output-json-file",
        output_json_file,
        "--output-path",
        tmp_path,
        "--sequencing-technology",
        sequencing_technology,
        "--notifications-path",
        tmp_path,
    ]

    if ncov_csv:
        args.extend(
            [
                "--ncov-qc-csv-file",
                test_data_path / "pipeline_results_files" / ncov_csv,
            ]
        )
    rv = CliRunner().invoke(
        generate_pipeline_results_files,
        args,
    )

    assert rv.exit_code == 0

    expected_output_df = pd.read_csv(test_data_path / "pipeline_results_files" / exp_results_csv)
    calculated_output_df = pd.read_csv(output_csv_file)

    calculated_cols = calculated_output_df.columns
    expected_cols = expected_output_df.columns
    assert sorted(calculated_cols) == sorted(expected_cols)

    calculated_output_df.sort_values(by=[SAMPLE_ID], inplace=True)
    calculated_output_df.sort_index(axis=1, inplace=True)
    expected_output_df.sort_values(by=[SAMPLE_ID], inplace=True)
    expected_output_df.sort_index(axis=1, inplace=True)

    assert_frame_equal(
        calculated_output_df.reset_index(drop=True),
        expected_output_df.reset_index(drop=True),
    )

    if ncov_csv:
        check_sample_list(Path(tmp_path / SAMPLES_UNKNOWN_NCOV_QC_FILE), exp_lists[UNKNOWN_NCOV])
        check_sample_list(Path(tmp_path / SAMPLES_FAILED_NCOV_QC_FILE), exp_lists[FAILED_NCOV])
        check_sample_list(Path(tmp_path / SAMPLES_PASSED_NCOV_QC_FILE), exp_lists[PASSED_NCOV])

    check_sample_list(Path(tmp_path / SAMPLES_UNKNOWN_PANGOLIN_FILE), exp_lists[UNKNOWN_PANGOLIN])
    check_sample_list(Path(tmp_path / SAMPLES_FAILED_PANGOLIN_FILE), exp_lists[FAILED_PANGOLIN])
    check_sample_list(Path(tmp_path / SAMPLES_PASSED_PANGOLIN_FILE), exp_lists[PASSED_PANGOLIN])

    # compare the resultfiles.json
    with open(test_data_path / "pipeline_results_files" / exp_resultfiles_json) as json_fd:
        exp_resultfiles_json_dict = json.load(json_fd)
    with open(output_json_file) as json_fd:
        calc_resultfiles_json_dict = json.load(json_fd)

    exp_result_files_json_dict_full_path = {
        sample_id: sorted([f"{str(tmp_path)}/{f}" for f in files])
        for sample_id, files in exp_resultfiles_json_dict.items()
    }
    calc_result_files_json_dict_sorted = {
        sample_id: sorted(files) for sample_id, files in calc_resultfiles_json_dict.items()
    }

    assert exp_result_files_json_dict_full_path == calc_result_files_json_dict_sorted
