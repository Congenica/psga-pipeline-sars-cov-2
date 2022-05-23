from pathlib import Path
from ftplib import FTP

import click


SUBMIT_COMPLETE_MARK_FILE_NAME = "submit.ready"


def upload_file(ftp_session: FTP, file_to_upload: Path, destination: str, binary: bool = False) -> None:
    """
    Upload a file to FTP
    """
    if binary:
        with open(file_to_upload, "rb") as f:
            ftp_session.storbinary(f"STOR {destination}", f)
    else:
        # Still need to open in "rb", as storlines expects a bytes-like object
        with open(file_to_upload, "rb") as f:
            ftp_session.storlines(f"STOR {destination}", f)


# pylint: disable=unused-argument
@click.command()
@click.option(
    "--input-xml",
    type=click.Path(file_okay=True, readable=True, exists=True),
    required=True,
    help="submission form to GenBank. File acts like an envelope to the submission and includes the "
    "necessary instructions on how to process this submission.",
)
@click.option(
    "--input-zip",
    type=click.Path(file_okay=True, readable=True, exists=True),
    required=True,
    help="A .zip archive file containing the genome data to be uploaded to the appropriate submission folder",
)
@click.option(
    "--submit-id",
    type=str,
    required=True,
    help="Unique identifier of the submission to GenBank",
)
@click.option(
    "--url",
    type=str,
    required=True,
    help="GenBank remote url",
)
@click.option(
    "--username",
    type=str,
    required=True,
    help="GenBank remote username",
)
@click.option(
    "--password",
    type=str,
    required=True,
    help="GenBank remote password",
)
@click.option(
    "--directory",
    type=str,
    required=True,
    help="GenBank storage directory. May be 'Test' or 'Production'",
)
def submit_genbank_files(
    input_xml: str, input_zip: str, submit_id: str, url: str, username: str, password: str, directory: str
) -> None:
    """
    Submit genome files to GenBank remote
    """
    signal_submit_complete_file = Path(SUBMIT_COMPLETE_MARK_FILE_NAME)
    signal_submit_complete_file.touch()

    session = FTP(host=f"{url}", user=username, passwd=password)
    session.cwd(directory)

    # Create and enter unique submission folder
    unique_upload_folder = submit_id
    if unique_upload_folder not in session.nlst():
        session.mkd(unique_upload_folder)
    session.cwd(unique_upload_folder)

    zip_file = Path(input_zip)
    upload_file(ftp_session=session, file_to_upload=zip_file, destination=zip_file.name, binary=True)
    print(f"Uploaded file {zip_file.name}")

    xml_file = Path(input_xml)
    upload_file(ftp_session=session, file_to_upload=xml_file, destination=xml_file.name)
    print(f"Uploaded file {xml_file.name}")

    upload_file(
        ftp_session=session, file_to_upload=signal_submit_complete_file, destination=signal_submit_complete_file.name
    )
    print(f"Uploaded final file {signal_submit_complete_file.name}. Upload is complete")

    session.quit()


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    submit_genbank_files()
