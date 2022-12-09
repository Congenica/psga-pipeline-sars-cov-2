from os.path import join as join_path  # used to join FS paths and S3 URIs
from pathlib import Path, PosixPath
from typing import Dict, List, Set
from functools import partial, reduce
import random
import csv
import json
from json import JSONEncoder

import click
import pandas as pd

from scripts.util.logger import get_structlog_logger
from scripts.util.metadata import EXPECTED_HEADERS as EXPECTED_METADATA_HEADERS, SAMPLE_ID, ILLUMINA, ONT
from scripts.validation.check_csv_columns import check_csv_columns

from collections import namedtuple

# These are all of the files that get put in the resultfiles.json
JsonOutFile = namedtuple('JsonOutFile', 'filepath typedescription')
OUTPUT_FILES = [
    # TODO where does GCF_004153365 come from?
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/logs/custom_dumpsoftwareversions/versions.yml', 'software-versions/customdump'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/software_versions.yml', 'software-versions/all'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/logs/assembly_qc/checkm.log', 'software-versions/checkm'),
    JsonOutFile('{output_path}/bactopia_one/{sample_id}_1.fastq.gz', 'input-reads-1'),
    JsonOutFile('{output_path}/bactopia_one/{sample_id}_2.fastq.gz', 'input-reads-2'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.raw.vcf.gz', 'variants/raw'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.annotated.vcf.gz', 'variants/annotated'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.filt.vcf.gz', 'variants/filtered'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.bam', 'variants/bam'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.subs.vcf.gz', 'variants/subs'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.vcf.gz', 'variants/vcf'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/GCF_004153365/{sample_id}.consensus.fa.gz', 'variants/consensus'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/mlst/default/blast/{sample_id}-blast.json', 'mlst/blast'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/antimicrobial-resistance/{sample_id}-protein-report.txt', 'antimicrobial-resistance/protein'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/antimicrobial-resistance/{sample_id}-gene-report.txt', 'antimicrobial-resistance/gene'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/assembly/checkm/bins/{sample_id}/hmmer.tree.txt', 'checkm/hmmer'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/assembly/checkm/bins/{sample_id}/genes.gff', 'checkm/genes'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/assembly/checkm/checkm-results.txt', 'checkm/results'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/bactopia/nf-reports/bactopia-trace.txt', 'bactopia-run-info'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{sample_id}/{sample_id}.bam.bai', 'variants/bam-bai'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{sample_id}/{sample_id}consensus.subs.fa.gz', 'variants/consensus-subs'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{sample_id}/{sample_id}consensus.subs.masked.fa.gz', 'variants/consensus.subs.masked'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/quality-control/summary/{sample_id}_R1-final_fastqc.zip', 'QC/R1-final-fastqc'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/quality-control/summary/{sample_id}_R2-final_fastqc.zip', 'QC/R2-final-fastqc'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/quality-control/summary/{sample_id}_R1-original_fastqc.zip', 'QC/R1-original-fastqc'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/quality-control/summary/{sample_id}_R2-original_fastqc.zip', 'QC/R2-original-fastqc'),
    JsonOutFile('{output_path}/bactopia_one/quast_assembly.zip', 'assembly/quast-zip'),
    #      TODO all local-single-sample/*/variants/*/*.vcf.gz
]


@click.command()
@click.option(
    "--metadata-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="the sample metadata file",
)
@click.option(
    "--output-csv-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file merging the results from ncov and pangolin",
)
@click.option(
    "--output-json-file",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="output file containing all the expected files per sample",
)
@click.option(
    "--output-path",
    type=str,
    required=True,
    help="output_path output path where sample result files are stored (e.g. s3://bucket/path/analysis_run)",
)
@click.option(
    "--sequencing-technology",
    type=click.Choice([ILLUMINA, ONT], case_sensitive=True),
    required=True,
    help="the sequencer technology used for sequencing the samples",
)
def generate_results(
    metadata_file: str,
    output_csv_file: str,
    output_json_file: str,
    output_path: str,
    sequencing_technology: str,
) -> None:
    """
    Generate pipeline results files
    """
    print('metadata_file: ' + metadata_file)
    print('output_csv_file: ' + output_csv_file)
    print('output_json_file: ' + output_json_file)
    print('output_path: ' + output_path)
    print('sequencing_technology: ' + sequencing_technology)

    # Make a list so we can iterate over this several times easily
    metadata = list(csv.DictReader(open(Path(metadata_file))))

    with open(output_csv_file, 'w', encoding='UTF8') as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow(['SAMPLE_ID', 'STATUS'])
        for sample in metadata:
            # For the demo everything passes TODO come back to this
            writer.writerow([sample['SAMPLE_ID'], 'PASSED'])

    json_out = dict()

    # Build resultfiles.json
    for sample in metadata:
        out_files = list()
        for of in OUTPUT_FILES:
            # Make sure file and type not already in list (psga framework required no duplicates)
            # TODO in production version error rather than silently fail
            file_in_list = False
            f = of.filepath.format(output_path=output_path, sample_id=sample['SAMPLE_ID'])
            t = of.typedescription
            for item in out_files:
                if item["file"] == f or item["type"] == t:
                    print('file: ', item["file"])
                    print('type: ', item["type"])
                    print('f: ', f)
                    print('t: ', t)
                    file_in_list = True
            if not file_in_list:
                out_files.append({
                    "file": f,
                    "type": t
                })
            else:
                print('WARNING: found duplicate file {t} : {f}'.format(t=t, f=f))
        json_out[sample['SAMPLE_ID']] = out_files

    with open(output_json_file, 'w') as f:
        json.dump(json_out, f)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_results()
