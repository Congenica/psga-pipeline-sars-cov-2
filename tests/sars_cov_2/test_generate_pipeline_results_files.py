from pathlib import Path
import pytest
import json
from click.testing import CliRunner
import pandas as pd
from pandas.testing import assert_frame_equal

from scripts.sars_cov_2.generate_pipeline_results_files import generate_pipeline_results_files
from scripts.util.metadata import SAMPLE_ID


@pytest.mark.parametrize(
    "use_s3_uri",
    [False, True],
)
@pytest.mark.parametrize(
    "metadata,contamination_removal_csv,primer_autodetection_csv,ncov_csv,pangolin_csv,"
    "analysis_run_name,sequencing_technology,exp_results_csv,exp_resultfiles_json",
    [
        (
            "metadata_illumina.csv",
            "contamination_removal.csv",
            "primer_autodetection.csv",
            "ncov_test.qc.csv",
            "all_lineages_report.csv",
            "just_a_name",
            "illumina",
            "results_illumina_ont.csv",
            "resultfiles_illumina.json",
        ),
        (
            "metadata_ont.csv",
            "contamination_removal.csv",
            "primer_autodetection.csv",
            "ncov_test.qc.csv",
            "all_lineages_report.csv",
            "just_a_name",
            "ont",
            "results_illumina_ont.csv",
            "resultfiles_ont.json",
        ),
        (
            "metadata_unknown.csv",
            None,
            None,
            None,
            "all_lineages_report.csv",
            "just_a_name",
            "unknown",
            "results_unknown.csv",
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
    ncov_csv,
    pangolin_csv,
    analysis_run_name,
    sequencing_technology,
    exp_results_csv,
    exp_resultfiles_json,
    use_s3_uri,
):

    output_csv_file = Path(tmp_path / "results.csv")
    output_json_file = Path(tmp_path / "resultfiles.json")

    output_path = "s3://bucket/analysis_run" if use_s3_uri else tmp_path

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
        output_path,
        "--sequencing-technology",
        sequencing_technology,
    ]

    if ncov_csv:
        args.extend(
            [
                "--contamination-removal-csv-file",
                test_data_path / "pipeline_results_files" / contamination_removal_csv,
                "--primer-autodetection-csv-file",
                test_data_path / "pipeline_results_files" / primer_autodetection_csv,
                "--ncov-qc-csv-file",
                test_data_path / "pipeline_results_files" / ncov_csv,
            ]
        )
    rv = CliRunner().invoke(
        generate_pipeline_results_files,
        args,
    )

    assert rv.exit_code == 0

    # assert results.csv
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

    # assert resultfiles.json
    with open(test_data_path / "pipeline_results_files" / exp_resultfiles_json) as json_fd:
        exp_resultfiles_json_dict = json.load(json_fd)
    with open(output_json_file) as json_fd:
        calc_resultfiles_json_dict = json.load(json_fd)

    exp_result_files_json_dict_full_path = {
        sample_id: sorted([f"{str(output_path)}/{f}" for f in files])
        for sample_id, files in exp_resultfiles_json_dict.items()
    }
    calc_result_files_json_dict_sorted = {
        sample_id: sorted(files) for sample_id, files in calc_resultfiles_json_dict.items()
    }

    assert exp_result_files_json_dict_full_path == calc_result_files_json_dict_sorted
