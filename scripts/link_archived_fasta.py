from os import symlink
from os.path import basename, isfile, islink, join
from pathlib import Path

import click

FASTA_FILE_EXTENSION = "fasta"


@click.command()
@click.option(
    "--source",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    default="",
    envvar="COVID_PIPELINE_FASTA_PATH",
    help="Source directory to search for .fasta files",
)
@click.option(
    "--destination",
    type=click.Path(exists=True, dir_okay=True, readable=True),
    required=True,
    help="Destination directory to store the links",
)
def link_archived_fasta(source: str, destination: str) -> None:
    """
    Create symlinks in the destination directory to the archived FASTA files if not already present
    """
    for source_path in Path(source).rglob(f"*.{FASTA_FILE_EXTENSION}"):
        click.echo(f"Processing file: {source_path}")
        dest_path = join(destination, basename(source_path))
        if not isfile(dest_path) and not islink(dest_path):
            symlink(source_path, dest_path)
            click.echo("... file not found in the destination dir")
            click.echo("... created symlink in the destination dir")


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    link_archived_fasta()
