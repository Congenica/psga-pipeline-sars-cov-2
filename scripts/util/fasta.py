from pathlib import Path

from Bio import SeqIO


FASTA_FILE_HANDLE = "fasta"


def parse_sequences(filename: Path) -> list:
    return list(SeqIO.parse(filename, FASTA_FILE_HANDLE))


def merge_fasta(sequences: list, output_file: Path):
    with open(output_file, "w") as concatenated_fasta:
        SeqIO.write(sequences, concatenated_fasta, FASTA_FILE_HANDLE)
