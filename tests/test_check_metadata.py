import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from scripts.check_metadata import check_metadata
from utils_tests import read_samples_from_file


@pytest.mark.parametrize(
    "metadata_file,analysis_run_name,analysis_run_columns,"
    "load_missing_samples,valid_samples,invalid_samples,exit_code,exception_msg",
    [
        (
            "good_metadata_illumina_fastq.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            True,
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_illumina_bam.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "bam",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            True,
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_fasta.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fasta",
                "workflow": "none",
                "pipeline_version": "1.0.0",
            },
            True,
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_medaka_fastq.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "medaka_artic",
                "pipeline_version": "1.0.0",
            },
            True,
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            [],
            0,
            None,
        ),
        (
            "good_metadata_fasta.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fasta",
                "workflow": "medaka_artic",
                "pipeline_version": "1.0.0",
            },
            True,
            [],
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            1,
            "Error: workflow medaka_artic does not support input file type fasta\n",
        ),
        (
            "good_metadata_illumina_fastq.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            False,
            [],
            ["37a36d1c-5985-4836-87b5-b36bac75d81b", "985347c5-ff6a-454c-ac34-bc353d05dd70"],
            1,
            "Invalid row for sample 37a36d1c-5985-4836-87b5-b36bac75d81b:\n"
            + "Sample 37a36d1c-5985-4836-87b5-b36bac75d81b not found in the database, but listed in pipeline metadata\n"
            + "Invalid row for sample 985347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "Sample 985347c5-ff6a-454c-ac34-bc353d05dd70 not found in the database, but listed in pipeline metadata\n"
            + "Error: Errors encountered for sample ids: 37a36d1c-5985-4836-87b5-b36bac75d81b, "
            + "985347c5-ff6a-454c-ac34-bc353d05dd70\n",
        ),
        (
            "good_metadata_medaka_fastq.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "bam",
                "workflow": "medaka_artic",
                "pipeline_version": "1.0.0",
            },
            False,
            [],
            [],
            1,
            "Error: workflow medaka_artic does not support input file type bam\n",
        ),
        (
            "good_metadata_illumina_fastq.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "unknown",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            False,
            [],
            [],
            2,
            "Error: Invalid value for '--input-file-type'",
        ),
        (
            "good_metadata_illumina_fastq.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "fake_workflow",
                "pipeline_version": "1.0.0",
            },
            False,
            [],
            [],
            2,
            "Error: Invalid value for '--workflow'",
        ),
        (
            "bad_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            False,
            [],
            [
                '""',
                "#()aadd",
                "185347c5-ff6a-454c-ac34-bc353d05dd70",
                "27a36d1c-5985-4836-87b5-b36bac75d81b",
                "385347c5-ff6a-454c-ac34-bc353d05dd70",
                "485347c5-ff6a-454c-ac34-bc353d05dd70",
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
            + "Invalid row for sample 385347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "md5_1 for 385347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
            + "Invalid row for sample 485347c5-ff6a-454c-ac34-bc353d05dd70:\n"
            + "md5_2 for 485347c5-ff6a-454c-ac34-bc353d05dd70 not available\n"
            + "Error: Errors encountered for sample ids: , #()aadd, 185347c5-ff6a-454c-ac34-bc353d05dd70, "
            + "27a36d1c-5985-4836-87b5-b36bac75d81b, 385347c5-ff6a-454c-ac34-bc353d05dd70, "
            + "485347c5-ff6a-454c-ac34-bc353d05dd70\n",
        ),
    ],
)
def test_check_metadata(
    db_session,
    tmp_path,
    test_data_path,
    metadata_file,
    analysis_run_name,
    analysis_run_columns,
    load_missing_samples,
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
        "--pipeline-version",
        analysis_run_columns["pipeline_version"],
        "--primer-scheme-name",
        analysis_run_columns["primer_scheme_name"],
        "--primer-scheme-version",
        analysis_run_columns["primer_scheme_version"],
        "--input-file-type",
        analysis_run_columns["input_file_type"],
        "--workflow",
        analysis_run_columns["workflow"],
        "--samples-with-valid-metadata-file",
        valid_samples_path,
        "--samples-with-invalid-metadata-file",
        invalid_samples_path,
    ]

    if load_missing_samples:
        cmd_config.append("--load-missing-samples")

    rv = CliRunner().invoke(
        check_metadata,
        cmd_config,
    )

    print(rv.output)
    print(rv.exception)
    assert rv.exit_code == exit_code
    samples = db_session.query(Sample).all()

    if exit_code == 2:
        assert exception_msg in str(rv.output)

    else:
        if exit_code == 0:
            assert len(samples) == 2

            analysis_run = (
                db_session.query(AnalysisRun)
                .filter(
                    AnalysisRun.analysis_run_name == analysis_run_name,
                )
                .one_or_none()
            )
            assert analysis_run is not None
            for col_name, col_val in analysis_run_columns.items():
                col_val_str = str(col_val)
                if col_name in ("input_file_type", "workflow"):
                    col_val_str = col_val_str.upper()
                # get column of sample from string, dynamically
                assert str(getattr(analysis_run, col_name)) == col_val_str

        elif exit_code == 1:
            assert len(samples) == 0
            assert rv.output == exception_msg

        if analysis_run_columns["workflow"] != "medaka_artic" and analysis_run_columns["input_file_type"] != "bam":
            # these files are not generated for medaka_artic + bam, as this combination is not currently supported

            # check lists of valid and invalid samples
            processed_valid_samples = read_samples_from_file(valid_samples_path)
            assert sorted(valid_samples) == sorted(processed_valid_samples)

            processed_invalid_samples = read_samples_from_file(invalid_samples_path)
            assert sorted(invalid_samples) == sorted(processed_invalid_samples)
