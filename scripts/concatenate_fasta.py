from pathlib import Path

import click
from Bio import SeqIO

from reheader_fasta import FASTA_FILE_HANDLE


@click.command()
@click.argument("source")
@click.option("--output", type=str, required=True, help="File to which to write all fasta files")
@click.option("--root-genome", type=str, required=True, help="Root genome to be prepended when concatenating")
def concatenate_fasta(source: str, output: str, root_genome: str) -> None:
    def parse_sequences(filename):
        return list(SeqIO.parse(filename, FASTA_FILE_HANDLE))

    records = parse_sequences(root_genome)
    for fasta_file in Path(source).rglob(f"*.{FASTA_FILE_HANDLE}"):
        records.extend(parse_sequences(fasta_file))

    with open(output, "w") as concatenated_fasta:
        SeqIO.write(records, concatenated_fasta, FASTA_FILE_HANDLE)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    concatenate_fasta()
