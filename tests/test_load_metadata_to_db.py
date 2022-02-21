import pytest
from click.testing import CliRunner

from scripts.db.models import AnalysisRun, Sample
from scripts.load_metadata_to_db import load_metadata


@pytest.mark.parametrize(
    "metadata_file,analysis_run_name,analysis_run_columns,exit_code,exception_msg",
    [
        (
            "good_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            0,
            None,
        ),
        (
            "good_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "bam",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            0,
            None,
        ),
        (
            "good_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "medaka_artic",
                "pipeline_version": "1.0.0",
            },
            0,
            None,
        ),
        (
            "good_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "bam",
                "workflow": "medaka_artic",
                "pipeline_version": "1.0.0",
            },
            1,
            "Error: medaka_artic workflow does not support input bam files\n",
        ),
        (
            "good_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "unknown",
                "workflow": "illumina_artic",
                "pipeline_version": "1.0.0",
            },
            2,
            "Error: Invalid value for '--input-file-type'",
        ),
        (
            "good_metadata.tsv",
            "just_a_name",
            {
                "primer_scheme_name": "nCoV-2019",
                "primer_scheme_version": "V3",
                "input_file_type": "fastq",
                "workflow": "fake_workflow",
                "pipeline_version": "1.0.0",
            },
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
            1,
            "Invalid row for sample ID HAM44444:\n"
            + 'ASSIGN DATE "32/02/2020" is not a valid date: day is out of range for month\n'
            + "Error: Errors encountered: HAM44444\n",
        ),
    ],
)
def test_load_metadata_to_db(
    db_session, test_data_path, metadata_file, analysis_run_name, analysis_run_columns, exit_code, exception_msg
):
    rv = CliRunner().invoke(
        load_metadata,
        [
            "--file",
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
        ],
    )

    assert rv.exit_code == exit_code
    samples = db_session.query(Sample).all()

    if exit_code == 0:
        assert len(samples) == 6

        analysis_run = (
            db_session.query(AnalysisRun)
            .filter(
                AnalysisRun.analysis_run_name == analysis_run_name,
            )
            .one_or_none()
        )
        assert analysis_run is not None
        for col_name, col_val in analysis_run_columns.items():
            # get column of sample from string, dynamically
            assert str(getattr(analysis_run, col_name)) == str(col_val)

    elif exit_code == 1:
        assert len(samples) == 0
        assert rv.output == exception_msg

    else:
        assert exception_msg in str(rv.output)
