import pytest
from click.testing import CliRunner

from scripts.sars_cov_2.check_metadata import check_metadata, validate_metadata
from utils_tests import read_samples_from_file


@pytest.mark.parametrize(
    "metadata_file,sequencing_technology,valid_samples,invalid_samples",
    [
        (
            "good_metadata_illumina.csv",
            "illumina",
            [
                "37a36d1c-5985-4836-87b5-b36bac75d81b",
                "57a36d1c-5985-4836-87b5-b36bac75d81b",
                "885347c5-ff6a-454c-ac34-bc353d05dd70",
                "985347c5-ff6a-454c-ac34-bc353d05dd70",
                "47a36d1c-5985-4836-87b5-b36bac75d81b",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            [],
        ),
        (
            "bad_metadata.csv",
            "illumina",
            [
                "385347c5-ff6a-454c-ac34-bc353d05dd70",
                "#()aadd",
            ],
            [
                "",
                "185347c5-ff6a-454c-ac34-bc353d05dd70",
                "186647c5-ff6a-454c-ac34-bc353d05dd70",
                "27a36d1c-5985-4836-87b5-b36bac75d81b",
                "286647c5-ff6a-454c-ac34-bc353d05dd70",
            ],
        ),
    ],
)
def test_validate_metadata(
    test_data_path,
    metadata_file,
    sequencing_technology,
    valid_samples,
    invalid_samples,
):
    metadata_path = test_data_path / metadata_file

    samples = validate_metadata(metadata_path, sequencing_technology)

    assert sorted(samples.valid) == sorted(valid_samples)
    assert sorted(samples.invalid) == sorted(invalid_samples)


@pytest.mark.parametrize(
    "metadata_file,analysis_run_name,sequencing_technology,valid_samples,invalid_samples,exit_code,exception_msg",
    [
        # CORRECT METADATA
        (
            "good_metadata_illumina.csv",
            "just_a_name",
            "illumina",
            [
                "37a36d1c-5985-4836-87b5-b36bac75d81b",
                "57a36d1c-5985-4836-87b5-b36bac75d81b",
                "885347c5-ff6a-454c-ac34-bc353d05dd70",
                "985347c5-ff6a-454c-ac34-bc353d05dd70",
                "47a36d1c-5985-4836-87b5-b36bac75d81b",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            [],
            0,
            None,
        ),
        (
            "good_metadata_ont.csv",
            "just_a_name",
            "ont",
            [
                "37a36d1c-5985-4836-87b5-b36bac75d81b",
                "57a36d1c-5985-4836-87b5-b36bac75d81b",
                "885347c5-ff6a-454c-ac34-bc353d05dd70",
                "985347c5-ff6a-454c-ac34-bc353d05dd70",
                "47a36d1c-5985-4836-87b5-b36bac75d81b",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            [],
            0,
            None,
        ),
        (
            "good_metadata_fasta.csv",
            "just_a_name",
            "unknown",
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        # CORRECT METADATA, INCORRECT sequencing_technology
        (
            "good_metadata_illumina.csv",
            "just_a_name",
            "fake_technology",
            [],
            [],
            2,
            "Error: Invalid value for '--sequencing-technology'",
        ),
        (
            "good_metadata_illumina.csv",
            "just_a_name",
            "ont",
            [
                "37a36d1c-5985-4836-87b5-b36bac75d81b",
                "47a36d1c-5985-4836-87b5-b36bac75d81b",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
                "57a36d1c-5985-4836-87b5-b36bac75d81b",
                "885347c5-ff6a-454c-ac34-bc353d05dd70",
                "985347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            [],
            0,
            None,
        ),
        (
            "good_metadata_ont.csv",
            "just_a_name",
            "illumina",
            ["47a36d1c-5985-4836-87b5-b36bac75d81b", "485347c5-ff6a-454c-ac34-bc353d05dd70"],
            [
                "37a36d1c-5985-4836-87b5-b36bac75d81b",
                "57a36d1c-5985-4836-87b5-b36bac75d81b",
                "885347c5-ff6a-454c-ac34-bc353d05dd70",
                "985347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            1,
            "Invalid row for SAMPLE_ID 37a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SEQ_FILE_2 for 37a36d1c-5985-4836-87b5-b36bac75d81b not available\n"
            + "Invalid row for SAMPLE_ID 57a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SEQ_FILE_2 for 57a36d1c-5985-4836-87b5-b36bac75d81b not available\n"
            + "Invalid row for SAMPLE_ID 885347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SEQ_FILE_2 for 885347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
            + "Invalid row for SAMPLE_ID 985347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SEQ_FILE_2 for 985347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
            + "Error: Errors encountered for sample ids: "
            + "37a36d1c-5985-4836-87b5-b36bac75d81b, 57a36d1c-5985-4836-87b5-b36bac75d81b, "
            + "885347c5-ff6a-454c-ac34-bc353d05dd70, 985347c5-ff6a-454c-ac34-bc353d05dd70\n",
        ),
        (
            "good_metadata_ont.csv",
            "just_a_name",
            "unknown",
            [],
            [
                "37a36d1c-5985-4836-87b5-b36bac75d81b",
                "47a36d1c-5985-4836-87b5-b36bac75d81b",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
                "57a36d1c-5985-4836-87b5-b36bac75d81b",
                "885347c5-ff6a-454c-ac34-bc353d05dd70",
                "985347c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            1,
            "Invalid row for SAMPLE_ID 37a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SAMPLE_ID: 37a36d1c-5985-4836-87b5-b36bac75d81b has invalid file for sequencing technology unknown. "
            + "Supported files are ['fasta', 'fa']\n"
            + "Invalid row for SAMPLE_ID 57a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SAMPLE_ID: 57a36d1c-5985-4836-87b5-b36bac75d81b has invalid file for sequencing technology unknown. "
            + "Supported files are ['fasta', 'fa']\n"
            + "Invalid row for SAMPLE_ID 885347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SAMPLE_ID: 885347c5-ff6a-454c-ac34-bc353d05dd70 has invalid file for sequencing technology unknown. "
            + "Supported files are ['fasta', 'fa']\n"
            + "Invalid row for SAMPLE_ID 985347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SAMPLE_ID: 985347c5-ff6a-454c-ac34-bc353d05dd70 has invalid file for sequencing technology unknown. "
            + "Supported files are ['fasta', 'fa']\n"
            + "Invalid row for SAMPLE_ID 47a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SAMPLE_ID: 47a36d1c-5985-4836-87b5-b36bac75d81b has invalid file for sequencing technology unknown. "
            + "Supported files are ['fasta', 'fa']\n"
            + "Invalid row for SAMPLE_ID 485347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SAMPLE_ID: 485347c5-ff6a-454c-ac34-bc353d05dd70 has invalid file for sequencing technology unknown. "
            + "Supported files are ['fasta', 'fa']\n"
            + "Error: Errors encountered for sample ids: "
            + "37a36d1c-5985-4836-87b5-b36bac75d81b, 47a36d1c-5985-4836-87b5-b36bac75d81b, "
            + "485347c5-ff6a-454c-ac34-bc353d05dd70, 57a36d1c-5985-4836-87b5-b36bac75d81b, "
            + "885347c5-ff6a-454c-ac34-bc353d05dd70, 985347c5-ff6a-454c-ac34-bc353d05dd70\n",
        ),
        (
            "good_metadata_fasta.csv",
            "just_a_name",
            "ont",
            [],
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            1,
            "Invalid row for SAMPLE_ID 37a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SAMPLE_ID: 37a36d1c-5985-4836-87b5-b36bac75d81b has invalid file for sequencing technology ont. "
            + "Supported files are ['fastq', 'fastq.gz', 'fq', 'fq.gz', 'bam']\n"
            + "Invalid row for SAMPLE_ID 985347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SAMPLE_ID: 985347c5-ff6a-454c-ac34-bc353d05dd70 has invalid file for sequencing technology ont. "
            + "Supported files are ['fastq', 'fastq.gz', 'fq', 'fq.gz', 'bam']\n"
            + "Error: Errors encountered for sample ids: "
            + "37a36d1c-5985-4836-87b5-b36bac75d81b, 985347c5-ff6a-454c-ac34-bc353d05dd70\n",
        ),
        # INCORRECT METADATA
        (
            "bad_metadata.csv",
            "invalid_rows",
            "illumina",
            [
                "385347c5-ff6a-454c-ac34-bc353d05dd70",
                "#()aadd",
            ],
            [
                '""',
                "185347c5-ff6a-454c-ac34-bc353d05dd70",
                "186647c5-ff6a-454c-ac34-bc353d05dd70",
                "27a36d1c-5985-4836-87b5-b36bac75d81b",
                "286647c5-ff6a-454c-ac34-bc353d05dd70",
            ],
            1,
            "Invalid row for SAMPLE_ID :\n"
            + "SAMPLE_ID not available\n"
            + "Invalid row for SAMPLE_ID 185347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SEQ_FILE_1 for 185347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
            + "Invalid row for SAMPLE_ID 186647c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SAMPLE_ID: 186647c5-ff6a-454c-ac34-bc353d05dd70 has invalid file for sequencing technology illumina. "
            + "Supported files are ['fastq', 'fastq.gz', 'fq', 'fq.gz', 'bam']\n"
            + "Invalid row for SAMPLE_ID 27a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "SEQ_FILE_2 for 27a36d1c-5985-4836-87b5-b36bac75d81b not available\n"
            + "Invalid row for SAMPLE_ID 286647c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "SEQ_FILE_1 and SEQ_FILE_2 for 286647c5-ff6a-454c-ac34-bc353d05dd70 have different file types\n"
            + "Error: Errors encountered for sample ids: , "
            + "185347c5-ff6a-454c-ac34-bc353d05dd70, 186647c5-ff6a-454c-ac34-bc353d05dd70, "
            + "27a36d1c-5985-4836-87b5-b36bac75d81b, 286647c5-ff6a-454c-ac34-bc353d05dd70\n",
        ),
    ],
)
def test_check_metadata(
    tmp_path,
    test_data_path,
    metadata_file,
    analysis_run_name,
    sequencing_technology,
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
        "--sequencing-technology",
        sequencing_technology,
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

        if valid_samples is not None:
            # check lists of valid and invalid samples
            processed_valid_samples = read_samples_from_file(valid_samples_path)
            assert sorted(valid_samples) == sorted(processed_valid_samples)

        if invalid_samples is not None:
            processed_invalid_samples = read_samples_from_file(invalid_samples_path)
            assert sorted(invalid_samples) == sorted(processed_invalid_samples)
