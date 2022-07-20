import pytest
from click.testing import CliRunner

from scripts.sars_cov_2.check_metadata import check_metadata
from utils_tests import read_samples_from_file


@pytest.mark.parametrize(
    "metadata_file,analysis_run_name,analysis_run_columns," "valid_samples,invalid_samples,exit_code,exception_msg",
    [
        (
            "good_metadata_illumina_fastq.csv",
            "just_a_name",
            {
                "input_file_type": "unknown",
                "ncov_workflow": "illumina_artic",
            },
            [],
            [],
            2,
            "Error: Invalid value for '--input-file-type'",
        ),
        (
            "good_metadata_illumina_fastq.csv",
            "just_a_name",
            {
                "input_file_type": "fastq",
                "ncov_workflow": "fake_workflow",
            },
            [],
            [],
            2,
            "Error: Invalid value for '--ncov-workflow'",
        ),
        (
            "good_metadata_illumina_fastq.csv",
            "just_a_name",
            {
                "input_file_type": "fastq",
                "ncov_workflow": "illumina_artic",
            },
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_illumina_bam.csv",
            "just_a_name",
            {
                "input_file_type": "bam",
                "ncov_workflow": "illumina_artic",
            },
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_medaka_fastq.csv",
            "just_a_name",
            {
                "input_file_type": "fastq",
                "ncov_workflow": "medaka_artic",
            },
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_fasta.csv",
            "just_a_name",
            {
                "input_file_type": "fasta",
                "ncov_workflow": "no_ncov",
            },
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_illumina_fastq.csv",
            "invalid_workflow_filetype",
            {
                "input_file_type": "fasta",
                "ncov_workflow": "illumina_artic",
            },
            [],
            [],
            1,
            "Error: ncov workflow 'illumina_artic' does not support input file type 'fasta'\n",
        ),
        (
            "good_metadata_medaka_fastq.csv",
            "invalid_workflow_filetype",
            {
                "input_file_type": "bam",
                "ncov_workflow": "medaka_artic",
            },
            [],
            [],
            1,
            "Error: ncov workflow 'medaka_artic' does not support input file type 'bam'\n",
        ),
        (
            "good_metadata_medaka_fastq.csv",
            "invalid_workflow_filetype",
            {
                "input_file_type": "fasta",
                "ncov_workflow": "medaka_artic",
            },
            [],
            [],
            1,
            "Error: ncov workflow 'medaka_artic' does not support input file type 'fasta'\n",
        ),
        (
            "good_metadata_fasta.csv",
            "invalid_workflow_filetype",
            {
                "input_file_type": "fastq",
                "ncov_workflow": "no_ncov",
            },
            [],
            [],
            1,
            "Error: ncov workflow 'no_ncov' does not support input file type 'fastq'\n",
        ),
        (
            "good_metadata_fasta.csv",
            "invalid_workflow-filetype",
            {
                "input_file_type": "bam",
                "ncov_workflow": "no_ncov",
            },
            [],
            [],
            1,
            "Error: ncov workflow 'no_ncov' does not support input file type 'bam'\n",
        ),
        (
            "bad_metadata.csv",
            "invalid_rows",
            {
                "input_file_type": "fastq",
                "ncov_workflow": "illumina_artic",
            },
            [],
            [
                '""',
                "#()aadd",
                "185347c5-ff6a-454c-ac34-bc353d05dd70",
                "27a36d1c-5985-4836-87b5-b36bac75d81b",
            ],
            1,
            "Invalid row for sample :\n"
            + "sample_id not available\n"
            + "Invalid row for sample #()aadd:\n"
            + 'sample_id "#()aadd" is not a UUID\n'
            + "Invalid row for sample 185347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "file_1 for 185347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
            + "Invalid row for sample 27a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "file_2 for 27a36d1c-5985-4836-87b5-b36bac75d81b not available\n"
            + "Error: Errors encountered for sample ids: , #()aadd, 185347c5-ff6a-454c-ac34-bc353d05dd70, "
            + "27a36d1c-5985-4836-87b5-b36bac75d81b\n",
        ),
    ],
)
def test_check_metadata(
    tmp_path,
    test_data_path,
    metadata_file,
    analysis_run_name,
    analysis_run_columns,
    valid_samples,
    invalid_samples,
    exit_code,
    exception_msg,
):

    valid_samples_path = tmp_path / "valid_samples.txt"
    invalid_samples_path = tmp_path / "invalid_samples.txt"

    cmd_config = [
        "--metadata-path",
        test_data_path / metadata_file,
        "--analysis-run-name",
        analysis_run_name,
        "--input-file-type",
        analysis_run_columns["input_file_type"],
        "--ncov-workflow",
        analysis_run_columns["ncov_workflow"],
        "--samples-with-valid-metadata-file",
        valid_samples_path,
        "--samples-with-invalid-metadata-file",
        invalid_samples_path,
    ]

    rv = CliRunner().invoke(
        check_metadata,
        cmd_config,
    )

    assert rv.exit_code == exit_code

    if exit_code == 2:
        assert exception_msg in str(rv.output)

    else:
        if exit_code == 1:
            assert rv.output == exception_msg

        if valid_samples:
            # check lists of valid and invalid samples
            processed_valid_samples = read_samples_from_file(valid_samples_path)
            assert sorted(valid_samples) == sorted(processed_valid_samples)

        if invalid_samples:
            processed_invalid_samples = read_samples_from_file(invalid_samples_path)
            assert sorted(invalid_samples) == sorted(processed_invalid_samples)
