import re
from os.path import join
from pathlib import Path

import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTA_FILE_EXTENSION = "fa"
FASTA_FILE_HANDLE = "fasta"
SEQUENCE_DESCRIPTION = ""  # PSGA. Don't add this as it causes Nextstrain to crash
ASSEMBLY_LENGTHS_FILENAME = "sequence_lengths.text"


def convert_file(source_file: Path, output_dir: Path) -> None:
    sequence_lengths = {}
    for record in SeqIO.parse(source_file, FASTA_FILE_HANDLE):

        # ncov-illumina workflow
        # >Consensus_ERR4157960.primertrimmed.consensus_threshold_0.75_quality_20
        # ncov-nanopore workflow
        # >20200311_1427_X1_FAK72834_a3787181_barcode07/ARTIC/medaka
        sample_name = re.search(r"^(Consensus_)?[\w-]+", record.id)

        if not sample_name:
            click.echo(f'sample name not found in {source_file} in header "{record.id}" skipping')
            continue

        # remove string if present (e.g. ncov-illumina workflow)
        sample_name = sample_name.group(0).replace("Consensus_", "", 1)

        sequence_lengths[sample_name] = len(record.seq)
        output_file = join(output_dir, f"{sample_name}.{FASTA_FILE_HANDLE}")
        with open(output_file, "w") as output:
            new_record = SeqRecord(record.seq, id=sample_name, description=SEQUENCE_DESCRIPTION)
            SeqIO.write(new_record, output, FASTA_FILE_HANDLE)

    # write gene assembly lengths to file, for use in Nextstrain
    with open(join(output_dir, ASSEMBLY_LENGTHS_FILENAME), "a") as lengths_file:
        for sequence, length in sequence_lengths.items():
            lengths_file.write(f"{sequence}\t{length}\n")


@click.command()
# source directory to search for .fa files
@click.argument("source")
# directory to output. If not specified, will output to source directory
@click.argument("destination", default="")
def reheader_fasta(source: str, destination: str) -> None:
    """
    Genome sequences produce by ncov have sequence identifiers that include
    QC parameters. This is to get rid of them.
    We only consider consensus fasta files
    """
    destination = destination or source

    # ncov returns consensus *.fa files when executed with the illumina workflow, but
    # it returns consensus *.fasta files when executed with the nanopore workflow
    # Here we rename any *.fasta to *.fa as 'fasta' is the extension of our output files
    for path in Path(source).rglob(f"*.consensus.{FASTA_FILE_HANDLE}"):
        path.rename(path.with_suffix(f".{FASTA_FILE_EXTENSION}"))

    for path in Path(source).rglob(f"*.consensus.{FASTA_FILE_EXTENSION}"):
        click.echo(f"processing file {path}")
        convert_file(path, Path(destination))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    reheader_fasta()
