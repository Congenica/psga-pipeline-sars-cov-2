from typing import List
from pathlib import Path
import csv
import zipfile

import click

from scripts.db.database import session_handler
from scripts.db.models import Sample
from scripts.util.fasta import FASTA_FILE_HANDLE, parse_sequences, merge_fasta
from scripts.genbank.submission import (
    Submission,
    typeOrganization,
    typeType,
    roleType,
    ActionType,
    DescriptionType,
    AddFilesType,
    typeTargetDb,
    FileType,
    DataTypeType,
    typeFileAttribute,
    typeIdentifier,
    typeSPUID,
    typeName,
)

GENBANK_WIZARD_VALUE = "BankIt_SARSCoV2_api"
SEQUENCE_FSA_FILE_NAME = "sequence.fsa"
SOURCE_METADATA_TABLE_FILE_NAME = "source.src"
SUBMISSION_TEMPLATE_FILE_NAME = "template.sbt"
SOURCE_METADATA_ORGANISM = "Severe acute respiratory syndrome coronavirus 2"
SOURCE_METADATA_ISOLATE_VIRUS = "SARS-CoV-2"
SOURCE_METADATA_ISOLATE_HOST = (
    "human"  # common or scientific name of the host animal from which the virus was located, for example: Homo sapiens
)
SOURCE_METADATA_ISOLATE_COUNTRY_ABBREVIATION = "BHR"  # virus host country three letter abbreviation
SOURCE_METADATA_HOST = "Homo sapiens"
SOURCE_METADATA_COUNTRY = "Bahrain"  # the country where the sample was isolated.
# See https://www.ncbi.nlm.nih.gov/genbank/collab/country/ for INSDC country list.
# Additional locality information can be included after the colon in the country, for example: USA: Maryland

SOURCE_METADATA_COLUMNS = ["sequence_ID", "organism", "isolate", "host", "collection-date", "country"]


def get_unsubmitted_sample_names() -> List[str]:
    """
    Return a list of all samples names available, which are not submitted to GenBank yet
    """
    with session_handler() as session:
        samples = session.query(Sample).filter(Sample.genbank_submit_id.is_(None)).all()
        return [sample.lab_id for sample in samples]


def merge_fastas_to_sequence_fsa(fasta_files: List[Path], output_file: Path) -> Path:
    """
    Merge all required fastas to single file
    """
    records = []
    for fasta_file in fasta_files:
        records.extend(parse_sequences(fasta_file))

    merge_fasta(sequences=records, output_file=output_file)
    return output_file


def create_metadata_table_file(sample_names: List[str], output_file: Path) -> Path:
    """
    Metadata .tsv file, describing information about sequences included in the submission
    """
    with session_handler() as session:
        samples = session.query(Sample).filter(Sample.lab_id.in_(sample_names)).all()
        metadata_table = [
            {
                "sequence_ID": sample.lab_id,
                "organism": SOURCE_METADATA_ORGANISM,
                "isolate": f"{SOURCE_METADATA_ISOLATE_VIRUS}/"
                f"{SOURCE_METADATA_ISOLATE_HOST}/"
                f"{SOURCE_METADATA_ISOLATE_COUNTRY_ABBREVIATION}/"
                f"{sample.lab_id}/"
                f"{sample.date_collected.strftime('%Y')}",
                "host": SOURCE_METADATA_HOST,
                "collection-date": sample.date_collected.strftime("%Y-%m-%d"),
                "country": SOURCE_METADATA_COUNTRY,
            }
            for sample in samples
        ]
    with open(output_file, "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        header_keys = SOURCE_METADATA_COLUMNS
        tsv_writer.writerow(header_keys)
        for sequence_record in metadata_table:
            tsv_writer.writerow([sequence_record[key] for key in header_keys])

    return output_file


def generate_submission(
    submit_name: str, submitter: str, zip_file_name: str, spuid_namespace: str, spuid_unique_value: str
) -> Submission:
    """
    Create submission xml, providing details about the author and type of submission
    """
    submission = Submission(
        Description=DescriptionType(
            Comment=submit_name,
            Organization=[typeOrganization(type_=typeType.CENTER, role=roleType.OWNER, Name=typeName(submitter))],
        ),
        Action=[
            ActionType(
                AddFiles=AddFilesType(
                    target_db=typeTargetDb.GEN_BANK,
                    File=[FileType(file_path=zip_file_name, DataType=DataTypeType.GENBANKSUBMISSIONPACKAGE)],
                    Attribute=[typeFileAttribute(name="wizard", valueOf_=GENBANK_WIZARD_VALUE)],
                    Identifier=typeIdentifier(
                        SPUID=typeSPUID(spuid_namespace=spuid_namespace, valueOf_=spuid_unique_value)
                    ),
                )
            )
        ],
    )
    return submission


def create_submission_archive_file(
    sample_names: List[str],
    sample_sequence_files: List[Path],
    template: Path,
    output_sequence_data: Path,
    output_source_metadata_table: Path,
    output_submission_archive: Path,
) -> Path:
    """
    Create one of two mandatory submission files - archive with all required sequence data and metadata supporting it
    """
    merged_fasta_file = merge_fastas_to_sequence_fsa(
        fasta_files=sample_sequence_files, output_file=output_sequence_data
    )
    source_metadata_table_file = create_metadata_table_file(
        sample_names=sample_names, output_file=output_source_metadata_table
    )
    with zipfile.ZipFile(output_submission_archive, "w") as output_zip:
        output_zip.write(merged_fasta_file, arcname=SEQUENCE_FSA_FILE_NAME)
        output_zip.write(source_metadata_table_file, arcname=SOURCE_METADATA_TABLE_FILE_NAME)
        output_zip.write(template, arcname=SUBMISSION_TEMPLATE_FILE_NAME)

    return output_submission_archive


def create_submission_metadata_file(
    submission_archive: Path,
    spuid_namespace: str,
    spuid_unique_value: str,
    submit_name: str,
    submitter: str,
    output_submission_metadata: Path,
) -> Path:
    """
    Create one of two mandatory submission files - metadata .xml file, describing information about submission and
    submitter
    """
    submission = generate_submission(
        submit_name=submit_name,
        submitter=submitter,
        zip_file_name=Path(submission_archive).name,
        spuid_namespace=spuid_namespace,
        spuid_unique_value=spuid_unique_value,
    )
    with open(output_submission_metadata, "w") as f:
        submission.export(f, 0, pretty_print=True)

    return output_submission_metadata


@click.command()
@click.option(
    "--input_sequence_fasta_directory",
    type=click.Path(dir_okay=True, readable=True),
    required=True,
    help="Directory, containing .fasta files. Sequence IDs must be unique for each sequence. "
    "Cannot contain spaces. may contain only the following characters - letters, digits, hyphens (-), "
    "underscores (_), periods (.), colons (:), asterisks (*), and number signs(#). should be under 25 characters."
    "For more information, see: https://submit.ncbi.nlm.nih.gov/genbank/help/#fasta",
)
@click.option(
    "--input_submission_template",
    type=click.Path(file_okay=True, readable=True, exists=True),
    required=True,
    help="(.sbt) Text file with submitter names and organizations, as well as publications associated with or "
    "describing the sequence. Users can generate submission template files at "
    "https://submit.ncbi.nlm.nih.gov/genbank/template/submission/. "
    "The saved template can be reused for multiple submissions.",
)
@click.option(
    "--output_sequence_data_fsa",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="Nucleotide sequences in FASTA format .fst",
)
@click.option(
    "--output_source_metadata_table_src",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="Tab-delimited text file .src, describing sequences within .fst file",
)
@click.option(
    "--output_submission_xml",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="Output file for writing GenBank submission xml",
)
@click.option(
    "--output_submission_zip",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="Output .zip file to put all genome information in",
)
@click.option(
    "--output_samples_submitted_file",
    type=click.Path(file_okay=True, writable=True),
    required=False,
    help="Output text file, which will be populated with a list of sample names for submission to GenBank",
)
@click.option("--submit_name", type=str, help="Name of submission")
@click.option("--submitter", type=str, help="Submitter name")
@click.option(
    "--spuid_namespace",
    type=str,
    help="Center/account abbreviation provided during account creation. "
    "This value remains the same for every submission",
)
@click.option(
    "--spuid_unique_value",
    type=str,
    help="Unique submission value. Must be unique for each submission from the submitter",
)
def generate_genbank_files(
    input_sequence_fasta_directory: str,
    input_submission_template: str,
    output_sequence_data_fsa: str,
    output_source_metadata_table_src: str,
    output_submission_xml: str,
    output_submission_zip: str,
    output_samples_submitted_file: str,
    submit_name: str,
    submitter: str,
    spuid_namespace: str,
    spuid_unique_value: str,
) -> None:
    """
    Create file collection, required for submitting genomes to GenBank
    """
    unsubmitted_sample_names = get_unsubmitted_sample_names()
    print(f"So far unsubmitted samples to GenBank: {', '.join(unsubmitted_sample_names)}")

    fasta_files = Path(input_sequence_fasta_directory).rglob(f"*.{FASTA_FILE_HANDLE}")
    fasta_names_and_files = {
        fasta_file.stem: fasta_file for fasta_file in fasta_files if fasta_file.stem in unsubmitted_sample_names
    }
    samples_to_submit = list(fasta_names_and_files.keys())
    samples_to_submit.sort()
    sample_files = [fasta_names_and_files[sample] for sample in samples_to_submit]
    print(f"Samples to submit: {', '.join(samples_to_submit)}")
    print(f"Samples files to submit: {', '.join(sample_file.name for sample_file in sample_files)}")

    submission_archive = create_submission_archive_file(
        sample_names=samples_to_submit,
        sample_sequence_files=sample_files,
        template=Path(input_submission_template),
        output_sequence_data=Path(output_sequence_data_fsa),
        output_source_metadata_table=Path(output_source_metadata_table_src),
        output_submission_archive=Path(output_submission_zip),
    )
    print(f"Done creating GenBank .zip file: {submission_archive.name}")

    submission_metadata = create_submission_metadata_file(
        submission_archive=submission_archive,
        spuid_namespace=spuid_namespace,
        spuid_unique_value=spuid_unique_value,
        submit_name=submit_name,
        submitter=submitter,
        output_submission_metadata=Path(output_submission_xml),
    )
    print(f"Done creating GenBank .xml file: {submission_metadata.name}")

    if output_samples_submitted_file:
        with open(output_samples_submitted_file, "w") as f:
            f.write("\n".join(samples_to_submit))
        print(f"Done writing samples being submitted to text file: {output_samples_submitted_file}")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_genbank_files()
