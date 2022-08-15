from pathlib import Path
import pytest
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.sars_cov_2.generate_pipeline_results_files import SAMPLE_ID_COL, generate_pipeline_results_files
from utils_tests import read_samples_from_file


def check_sample_list(input_path, expected_samples):
    processed_samples = read_samples_from_file(input_path)
    assert sorted(expected_samples) == sorted(processed_samples)


@pytest.mark.parametrize(
    "metadata,ncov_csv,pangolin_csv,analysis_run_name,exp_lists,exp_results_csv",
    [
        (
            "metadata.csv",
            "ncov_test.qc.csv",
            "all_lineages_report.csv",
            "just_a_name",
            {
                "ncov_unknown": ["985347c5-ff6a-454c-ac34-bc353d05dd70"],
                "ncov_failed": ["0774181d-fb20-4a73-b887-38af7eda9b38"],
                "ncov_passed": [
                    "e80f3c63-d139-4c14-bd72-7a43897ab40d",
                    "56a63f60-764d-4fc7-8764-a46023cbe324",
                    "a0951432-cd94-45b5-96d5-b721c037a451",
                ],
                "pangolin_unknown": ["a0951432-cd94-45b5-96d5-b721c037a451"],
                "pangolin_failed": ["e80f3c63-d139-4c14-bd72-7a43897ab40d"],
                "pangolin_passed": ["56a63f60-764d-4fc7-8764-a46023cbe324"],
            },
            "expected_merged_ncov_plus_pangolin.csv",
        ),
        (
            "metadata.csv",
            None,
            "all_lineages_report.csv",
            "just_a_name",
            {
                "pangolin_unknown": [
                    "985347c5-ff6a-454c-ac34-bc353d05dd70",
                    "0774181d-fb20-4a73-b887-38af7eda9b38",
                    "a0951432-cd94-45b5-96d5-b721c037a451",
                ],
                "pangolin_failed": ["e80f3c63-d139-4c14-bd72-7a43897ab40d"],
                "pangolin_passed": ["56a63f60-764d-4fc7-8764-a46023cbe324"],
            },
            "expected_merged_pangolin_only.csv",
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
    exp_lists,
    exp_results_csv,
):

    output_csv_file = Path(tmp_path / "results.csv")

    args = [
        "--analysis-run-name",
        analysis_run_name,
        "--metadata-file",
        test_data_path / "pipeline_results_files" / metadata,
        "--pangolin-csv-file",
        test_data_path / "pipeline_results_files" / pangolin_csv,
        "--output-csv-file",
        output_csv_file,
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

    calculated_output_df.sort_values(by=[SAMPLE_ID_COL], inplace=True)
    calculated_output_df.sort_index(axis=1, inplace=True)
    expected_output_df.sort_values(by=[SAMPLE_ID_COL], inplace=True)
    expected_output_df.sort_index(axis=1, inplace=True)

    assert_frame_equal(
        calculated_output_df.reset_index(drop=True),
        expected_output_df.reset_index(drop=True),
    )

    if ncov_csv:
        check_sample_list(Path(tmp_path / "samples_unknown_ncov_qc.txt"), exp_lists["ncov_unknown"])
        check_sample_list(Path(tmp_path / "samples_failed_ncov_qc.txt"), exp_lists["ncov_failed"])
        check_sample_list(Path(tmp_path / "samples_passed_ncov_qc.txt"), exp_lists["ncov_passed"])

    check_sample_list(Path(tmp_path / "samples_unknown_pangolin.txt"), exp_lists["pangolin_unknown"])
    check_sample_list(Path(tmp_path / "samples_failed_pangolin.txt"), exp_lists["pangolin_failed"])
    check_sample_list(Path(tmp_path / "samples_passed_pangolin.txt"), exp_lists["pangolin_passed"])
