
from os.path import join
from pathlib import Path

import click
from Bio import SeqIO

from reheader_fasta import FASTQ_FILE_HANDLE, SEQUENCE_DESCRIPTION

CONCATENATED_FASTA_FILENAME = f"nextstrain.{FASTQ_FILE_HANDLE}"


@click.command()
@click.argument("source", envvar="GENOME_FASTA_PATH")
def concatenate_fasta(source: str) -> None:
    def parse_sequences(filename):
        return [record for record in SeqIO.parse(filename, FASTQ_FILE_HANDLE)]

    # root repo dir
    root_genome_file = Path(__file__).parent.parent
    root_genome_file = root_genome_file / "data" / "FASTA" / f"{SEQUENCE_DESCRIPTION}.{FASTQ_FILE_HANDLE}"

    records = parse_sequences(root_genome_file)
    for fasta_file in Path(source).rglob(f"*.{FASTQ_FILE_HANDLE}"):
        records.extend(parse_sequences(fasta_file))

    with open(join(source, CONCATENATED_FASTA_FILENAME), "w") as concatenated_fasta:
        SeqIO.write(records, concatenated_fasta, FASTQ_FILE_HANDLE)


if __name__ == "__main__":
    concatenate_fasta()
