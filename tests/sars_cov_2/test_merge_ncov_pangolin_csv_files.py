from pathlib import Path
import pytest
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.sars_cov_2.merge_ncov_pangolin_csv_files import merge_ncov_pangolin_csv_files
from utils_tests import read_samples_from_file


def check_sample_list(input_path, expected_samples):
    processed_samples = read_samples_from_file(input_path)
    assert sorted(expected_samples) == sorted(processed_samples)


@pytest.mark.parametrize(
    "metadata,ncov_csv,pangolin_csv,analysis_run_name,exp_lists,exp_merged_output",
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
                "pangolin_unknown": ["0774181d-fb20-4a73-b887-38af7eda9b38", "a0951432-cd94-45b5-96d5-b721c037a451"],
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
def test_merge_ncov_pangolin_csv_files(
    tmp_path,
    test_data_path,
    metadata,
    ncov_csv,
    pangolin_csv,
    analysis_run_name,
    exp_lists,
    exp_merged_output,
):

    merged_output = Path(tmp_path / "output.csv")
    samples_unknown_ncov_qc = Path(tmp_path / "unknown_ncov.txt")
    samples_failed_ncov_qc = Path(tmp_path / "failed_ncov.txt")
    samples_passed_ncov_qc = Path(tmp_path / "passed_ncov.txt")
    samples_unknown_pangolin = Path(tmp_path / "unknown_pangolin.txt")
    samples_failed_pangolin = Path(tmp_path / "failed_pangolin.txt")
    samples_passed_pangolin = Path(tmp_path / "passed_pangolin.txt")

    args = [
        "--analysis-run-name",
        analysis_run_name,
        "--metadata-file",
        test_data_path / "merge_ncov_pangolin" / metadata,
        "--pangolin-csv-file",
        test_data_path / "merge_ncov_pangolin" / pangolin_csv,
        "--merged-output-csv-file",
        merged_output,
        "--samples-unknown-ncov-qc-file",
        samples_unknown_ncov_qc,
        "--samples-failed-ncov-qc-file",
        samples_failed_ncov_qc,
        "--samples-passed-ncov-qc-file",
        samples_passed_ncov_qc,
        "--samples-unknown-pangolin-file",
        samples_unknown_pangolin,
        "--samples-failed-pangolin-file",
        samples_failed_pangolin,
        "--samples-passed-pangolin-file",
        samples_passed_pangolin,
    ]

    if ncov_csv:
        args.extend(
            [
                "--ncov-qc-csv-file",
                test_data_path / "merge_ncov_pangolin" / ncov_csv,
            ]
        )
    rv = CliRunner().invoke(
        merge_ncov_pangolin_csv_files,
        args,
    )

    assert rv.exit_code == 0

    expected_output_df = pd.read_csv(test_data_path / "merge_ncov_pangolin" / exp_merged_output)
    calculated_output_df = pd.read_csv(merged_output)

    calculated_cols = calculated_output_df.columns
    expected_cols = expected_output_df.columns
    assert sorted(calculated_cols) == sorted(expected_cols)

    calculated_output_df.sort_values(by=["sample_id"], inplace=True)
    calculated_output_df.sort_index(axis=1, inplace=True)
    expected_output_df.sort_values(by=["sample_id"], inplace=True)
    expected_output_df.sort_index(axis=1, inplace=True)

    assert_frame_equal(
        calculated_output_df.reset_index(drop=True),
        expected_output_df.reset_index(drop=True),
    )

    if ncov_csv:
        check_sample_list(samples_unknown_ncov_qc, exp_lists["ncov_unknown"])
        check_sample_list(samples_failed_ncov_qc, exp_lists["ncov_failed"])
        check_sample_list(samples_passed_ncov_qc, exp_lists["ncov_passed"])

    check_sample_list(samples_unknown_pangolin, exp_lists["pangolin_unknown"])
    check_sample_list(samples_failed_pangolin, exp_lists["pangolin_failed"])
    check_sample_list(samples_passed_pangolin, exp_lists["pangolin_passed"])
