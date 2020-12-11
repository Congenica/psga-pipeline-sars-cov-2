import click


# pylint: disable=unused-argument
@click.command()
@click.option(
    "--input_xml",
    type=click.Path(file_okay=True, readable=True, exists=True),
    required=True,
    help="submission form to GenBank. File acts like an envelope to the submission and includes the "
    "necessary instructions on how to process this submission.",
)
@click.option(
    "--input_zip",
    type=click.Path(file_okay=True, readable=True, exists=True),
    required=True,
    help="A .zip archive file containing the genome data to be uploaded to the appropriate submission folder",
)
@click.option(
    "--submit_id",
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
def submit_genbank_files(
    input_xml: str, input_zip: str, submit_id: str, url: str, username: str, password: str
) -> None:
    """
    Submit genome files to GenBank remote
    """
    print("Not implemented yet!")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    submit_genbank_files()
