# genome sequences produce by ncov have sequence identifiers that include
# QC parameters. This is to get rid of them.

import os
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FASTQ_FILE_EXTENSION = "fa"
FASTQ_FILE_HANDLE = "fasta"
SEQUENCE_DESCRIPTION = "SARS-CoV-2"


def convert_file(source_file, output_dir):
    for record in SeqIO.parse(source_file, FASTQ_FILE_HANDLE):
        # >Consensus_ERR4157960.primertrimmed.consensus_threshold_0.75_quality_20
        sample_name = re.search(r"^Consensus_(\w+)", record.id)
        if not sample_name:
            print(
                f"sample name not found in {source_file} in header \"{record.id}\" skipping"
            )
            continue

        sample_name = sample_name.group(1)
        output_file = os.path.join(output_dir, f"{sample_name}.{FASTQ_FILE_HANDLE}")
        with open(output_file, "w") as output:
            new_record = SeqRecord(
                record.seq, id=sample_name, description=SEQUENCE_DESCRIPTION
            )
            SeqIO.write(new_record, output, FASTQ_FILE_HANDLE)


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
        convert_file(path, destination)
