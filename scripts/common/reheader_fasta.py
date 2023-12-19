from pathlib import Path

import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTA_FILE_EXTENSION = "fa"
FASTA_FILE_HANDLE = "fasta"
SEQUENCE_DESCRIPTION = ""  # PSGA. Don't add this as it causes Nextstrain to crash
ASSEMBLY_LENGTHS_FILENAME = "sequence_lengths.text"


def convert_file(source_file: Path, output_dir: Path) -> None:
    sequence_lengths = []
    # get the sample name from the filename
    source_file_extensions = "".join(source_file.suffixes)
    # use source_file.stem?
    sample_name = str(source_file.name).replace(source_file_extensions, "").rstrip("_1")

    output_file = output_dir / f"{sample_name}.{FASTA_FILE_HANDLE}"

    with open(output_file, "w") as output:
        for record in SeqIO.parse(source_file, FASTA_FILE_HANDLE):
            sequence_lengths.append(len(record.seq))
            new_record = SeqRecord(record.seq, id=sample_name, description=SEQUENCE_DESCRIPTION)
            SeqIO.write(new_record, output, FASTA_FILE_HANDLE)

    # write gene assembly lengths to file, for use in Nextstrain
    with open(output_dir / ASSEMBLY_LENGTHS_FILENAME, "a") as lengths_file:
        for sequence, length in enumerate(sequence_lengths):
            lengths_file.write(f"{sequence}\t{length}\n")


@click.command()
@click.option(
    "--input-dir",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Source directory to search for *.consensus.fa files",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(file_okay=True, writable=True),
    help="Output directory to write the reheadered fasta files",
)
def reheader_fasta(input_dir: str, output_dir: str) -> None:
    """
    Genome sequences produce by ncov have sequence identifiers that include
    QC parameters. This is to get rid of them.
    We only consider consensus fasta files
    """
    for path in Path(input_dir).rglob(f"*.consensus.{FASTA_FILE_EXTENSION}"):
        click.echo(f"processing file {path}")
        convert_file(path, Path(output_dir))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    reheader_fasta()
