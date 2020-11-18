# genome sequences produce by ncov have sequence identifiers that include
# QC parameters. This is to get rid of them.

import re
import sys
from os.path import join
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTQ_FILE_EXTENSION = "fa"
FASTQ_FILE_HANDLE = "fasta"
SEQUENCE_DESCRIPTION = "SARS-CoV-2"
ASSEMBLY_LENGTHS_FILENAME = "sequence_lengths.text"


def convert_file(source_file: Path, output_dir: Path):
    sequence_lengths = {}
    for record in SeqIO.parse(source_file, FASTQ_FILE_HANDLE):
        # >Consensus_ERR4157960.primertrimmed.consensus_threshold_0.75_quality_20
        sample_name = re.search(r"^Consensus_(\w+)", record.id)
        if not sample_name:
            print(
                f"sample name not found in {source_file} in header \"{record.id}\" skipping"
            )
            continue

        sample_name = sample_name.group(1)

        sequence_lengths[sample_name] = len(record.seq)
        output_file = join(output_dir, f"{sample_name}.{FASTQ_FILE_HANDLE}")
        with open(output_file, "w") as output:
            new_record = SeqRecord(
                record.seq, id=sample_name, description=SEQUENCE_DESCRIPTION
            )
            SeqIO.write(new_record, output, FASTQ_FILE_HANDLE)

    # write gene assembly lengths to file, for use in Nextstrain
    lengths_file = join(output_dir, ASSEMBLY_LENGTHS_FILENAME)
    with open(lengths_file, 'a') as lengths_file:
        for sequence, length in sequence_lengths.items():
            lengths_file.write(f'{sequence}\t{length}\n')


if __name__ == "__main__":
    # first parameter source dir, second (optional) destination dir
    if len(sys.argv) not in (2, 3):
        sys.exit(
            (
                f"Usage: {sys.argv[0]} <source directory to search "
                f"for .{FASTQ_FILE_EXTENSION} files> <directory to output>"
            )
        )

    destination = sys.argv[2] if len(sys.argv) == 3 else sys.argv[1]

    for path in Path(sys.argv[1]).rglob(f"*.{FASTQ_FILE_EXTENSION}"):
        print(f"processing file {path}")
        convert_file(path, Path(destination))
