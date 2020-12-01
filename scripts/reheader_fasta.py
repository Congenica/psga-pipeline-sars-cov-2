import re
from os.path import join
from pathlib import Path

import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTA_FILE_EXTENSION = "fa"
FASTA_FILE_HANDLE = "fasta"
SEQUENCE_DESCRIPTION = "SARS-CoV-2"
ASSEMBLY_LENGTHS_FILENAME = "sequence_lengths.text"


def convert_file(source_file: Path, output_dir: Path) -> None:
    sequence_lengths = {}
    for record in SeqIO.parse(source_file, FASTA_FILE_HANDLE):
        # >Consensus_ERR4157960.primertrimmed.consensus_threshold_0.75_quality_20
        sample_name = re.search(r"^Consensus_(\w+)", record.id)
        if not sample_name:
            click.echo(f'sample name not found in {source_file} in header "{record.id}" skipping')
            continue

        sample_name = sample_name.group(1)

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
@click.argument("destination", default="", envvar="COVID_PIPELINE_FASTA_PATH")
def reheader_fasta(source: str, destination: str) -> None:
    """Genome sequences produce by ncov have sequence identifiers that include
    QC parameters. This is to get rid of them.
    """
    destination = destination or source

    for path in Path(source).rglob(f"*.{FASTA_FILE_EXTENSION}"):
        click.echo(f"processing file {path}")
        convert_file(path, Path(destination))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    reheader_fasta()
