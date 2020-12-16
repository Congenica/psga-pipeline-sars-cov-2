from pathlib import Path

import click

from reheader_fasta import FASTA_FILE_HANDLE
from scripts.util.fasta import parse_sequences, merge_fasta


@click.command()
@click.argument("source")
@click.option("--output", type=str, required=True, help="File to which to write all fasta files")
@click.option("--root-genome", type=str, required=True, help="Root genome to be prepended when concatenating")
def concatenate_fasta(source: str, output: str, root_genome: str) -> None:
    records = parse_sequences(Path(root_genome))
    for fasta_file in Path(source).rglob(f"*.{FASTA_FILE_HANDLE}"):
        records.extend(parse_sequences(fasta_file))

    merge_fasta(records, Path(output))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    concatenate_fasta()
