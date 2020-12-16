import tempfile
from datetime import datetime
from zipfile import ZipFile

import pytest
from click.testing import CliRunner


from utils_tests import assert_files_are_equal
from scripts.db.models import Sample
from scripts.generate_genbank_files import generate_genbank_files


@pytest.fixture
def populated_db_session_with_samples(db_session):
    sample_names_not_submitted_yet = ["foo", "bar"]
    for sample in sample_names_not_submitted_yet:
        db_session.add(Sample(lab_id=sample, date_collected=datetime.fromtimestamp(0), metadata_loaded=True))

    sample_names_submitted = ["buzz"]
    for sample in sample_names_submitted:
        db_session.add(
            Sample(
                lab_id=sample, genbank_submit_id="foo", date_collected=datetime.fromtimestamp(0), metadata_loaded=True
            )
        )

    db_session.commit()
    yield db_session


@pytest.mark.parametrize(
    "fasta_dir,template,name,submitter,namespace,spuid,ref_fsa,ref_src,ref_txt,ref_xml",
    [
        (
            tempfile.TemporaryDirectory().name,
            "genbank.template.sbt",
            "Submission",
            "Submitter",
            "company",
            "id1",
            "empty.fsa",
            "empty.src",
            "empty.txt",
            "submission.xml",
        ),
        (
            "test_fasta",
            "genbank.template.sbt",
            "Submission",
            "Submitter",
            "company",
            "id1",
            "foo_bar.fsa",
            "foo_bar.src",
            "foo_bar.txt",
            "submission.xml",
        ),
    ],
)
def test_genbank_generating(
    populated_db_session_with_samples,
    test_data_path_genbank_input,
    test_data_path_genbank_reference,
    fasta_dir,
    template,
    name,
    submitter,
    namespace,
    spuid,
    ref_fsa,
    ref_src,
    ref_txt,
    ref_xml,
    tmp_path_factory,
):
    temp_dir = tmp_path_factory.mktemp("test_genbank_generating")
    input_dir = test_data_path_genbank_input
    ref_dir = test_data_path_genbank_reference
    out_fasta = temp_dir / "merged.fsa"
    out_metadata = temp_dir / "metadata.src"
    out_xml = temp_dir / "submission.xml"
    out_zip = temp_dir / "submission.zip"
    out_samples = temp_dir / "samples.txt"
    extract_dir = temp_dir / "extract"
    extract_dir.mkdir()

    rv = CliRunner().invoke(
        generate_genbank_files,
        [
            "--input_sequence_fasta_directory",
            input_dir / fasta_dir,
            "--input_submission_template",
            input_dir / template,
            "--output_sequence_data_fsa",
            out_fasta,
            "--output_source_metadata_table_src",
            out_metadata,
            "--output_submission_xml",
            out_xml,
            "--output_submission_zip",
            out_zip,
            "--output_samples_submitted_file",
            out_samples,
            "--submit_name",
            name,
            "--submitter",
            submitter,
            "--spuid_namespace",
            namespace,
            "--spuid_unique_value",
            spuid,
        ],
    )

    print(rv.output)
    assert rv.exit_code == 0

    assert_files_are_equal(out_xml, ref_dir / ref_xml)
    assert_files_are_equal(out_samples, ref_dir / ref_txt)
    with ZipFile(out_zip, "r") as submission_zip:
        submission_zip.extractall(extract_dir)

    assert_files_are_equal(extract_dir / "sequence.fsa", ref_dir / ref_fsa)
    assert_files_are_equal(extract_dir / "source.src", ref_dir / ref_src)
    assert_files_are_equal(extract_dir / "template.sbt", input_dir / template)
