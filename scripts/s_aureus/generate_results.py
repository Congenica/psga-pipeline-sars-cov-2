# flake8: noqa
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
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/logs/custom_dumpsoftwareversions/versions.yml', 'software-versions/customdump'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/software_versions.yml', 'software-versions/all'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/logs/assembly_qc/checkm.log', 'software-versions/checkm'),
    JsonOutFile('{output_path}/bactopia_one/{sample_id}_1.fastq.gz', 'input-reads-1'),
    JsonOutFile('{output_path}/bactopia_one/{sample_id}_2.fastq.gz', 'input-reads-2'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.raw.vcf.gz', 'variants/raw'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.annotated.vcf.gz', 'variants/annotated'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.filt.vcf.gz', 'variants/filtered'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.bam', 'variants/bam'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.subs.vcf.gz', 'variants/subs'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.vcf.gz', 'variants/vcf'),
    JsonOutFile('{output_path}/bactopia_one/local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.consensus.fa.gz', 'variants/consensus'),
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
]


def build_csv_list(sample_id, genome_assembly_name):
    """
    Gather all info for the CSV file
    Make a defined list of files that we are expecting.
    """
    CsvItem = namedtuple('CsvItem', 'header value')
    csv_list = list()
    # Add items to the list
    csv_list.append(CsvItem('SAMPLE_ID', sample_id))  # the key SAMPLE_ID must be in the results.csv
    csv_list.append(CsvItem('STATUS', 'Completed'))  # For the demo everything passes
    csv_list.append(CsvItem('QC_STATUS', 'pass'))  # The key QC_STATUS gets shown in the UI

    with open('local-single-sample/bactopia-tools/mykrobe/mykrobe/{sample_id}/logs/mykrobe/versions.yml'.format(sample_id=sample_id)) as software_versions:
        versions = list(csv.reader(software_versions, delimiter=":"))
        for row in versions:
            if 'checkm' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))
            if 'bactopia' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))
            if 'fastqc' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))

    with open('local-single-sample/bactopia-tools/mlst/mlst/software_versions.yml') as mlst_versions:
        m_versions = list(csv.reader(mlst_versions, delimiter=":"))
        for row in m_versions:
            if 'mlst' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))

    with open('local-single-sample/bactopia-tools/mykrobe/mykrobe/{sample_id}/logs/mykrobe/versions.yml'.format(sample_id=sample_id)) as mykrobe_versions:
        myk_versions = list(csv.reader(mykrobe_versions, delimiter=":"))
        for row in myk_versions:
            if 'mykrobe' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))

    with open('sample-{sample_id}-csv-files/annotation-summary.txt'.format(sample_id=sample_id)) as annotation_summary:
        annotations = csv.reader(annotation_summary, delimiter =":")
        for row in annotations:
            csv_list.append(CsvItem('Annotation-{}'.format(row[0]), row[1].lstrip().rstrip()))

    # Variant summary
    with open('local-single-sample/{sample_id}/variants/{genome_assembly_name}/{sample_id}.txt'.format(sample_id=sample_id, genome_assembly_name=genome_assembly_name)) as variant_summary_file:
        variant_summary = csv.reader(variant_summary_file, delimiter="\t")
        for row in variant_summary:
            if 'Variant' in row[0]:
                csv_list.append(CsvItem(row[0].rstrip().lstrip(), row[1].lstrip().rstrip()))

    with open('local-single-sample/{sample_id}/assembly/checkm/checkm-results.txt'.format(sample_id=sample_id)) as checkm_results_file:
        checkm_results = next(csv.DictReader(checkm_results_file, delimiter="\t"))
        csv_list.append(CsvItem('Marker lineage', checkm_results['Marker lineage']))
        csv_list.append(CsvItem('Completeness', checkm_results['Completeness']))
        csv_list.append(CsvItem('Contamination', checkm_results['Contamination']))
        csv_list.append(CsvItem('Strain heterogeneity', checkm_results['Strain heterogeneity']))

    with open('local-single-sample/{sample_id}/assembly/{sample_id}.json'.format(sample_id=sample_id)) as assembly_results_file:
        assembly_results = json.loads(assembly_results_file.read())
        for key, value in assembly_results.items():
            if key == 'fasta':
                continue
            csv_list.append(CsvItem(key , value))

    with open('local-single-sample/{sample_id}/antimicrobial-resistance/{sample_id}-protein-report.txt'.format(sample_id=sample_id)) as antimicrobial_protein_file:
        resistance_genes = list()
        for row in csv.DictReader(antimicrobial_protein_file, delimiter="\t"):
            resistance_genes.append(row['Gene symbol'])
        csv_list.append(CsvItem('Antimicrobial resistance genes', ' '.join(resistance_genes)))

    return csv_list



@click.command()
@click.option(
    "--output-path",
    type=str,
    required=True,
    help="output_path output path where sample result files are stored (e.g. s3://bucket/path/analysis_run)",
)
def generate_results(
    output_path: str,
) -> None:
    """
    Generate pipeline results files
    """
    # Get sample_id and genome_assembly_name
    with open('genome_assembly_name.txt') as genome_assembly_name_file:
        gaf = json.load(genome_assembly_name_file)
        genome_assembly_name = gaf['genome_assembly_name']
        sample_id = gaf['sample_id']
    print('genome_assembly_name: ' + genome_assembly_name)
    print('sample_id: ' + sample_id)

    csv_list = build_csv_list(sample_id, genome_assembly_name)

    with open('{sample_id}-results.csv'.format(sample_id=sample_id), 'w', encoding='UTF8') as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow([x.header for x in csv_list])
        writer.writerow([x.value for x in csv_list])

    # Build resultfiles.json
    json_out = dict()
    out_files = list()
    for of in OUTPUT_FILES:
        # Make sure file and type not already in list (psga framework required no duplicates)
        # TODO in production version error rather than fixing the problem with a warning
        file_in_list = False
        f = of.filepath.format(
            output_path=output_path, sample_id=sample_id, genome_assembly_name=genome_assembly_name)
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
    json_out[sample_id] = out_files

    with open('{sample_id}-resultfiles.json'.format(sample_id=sample_id), 'w') as f:
        json.dump(json_out, f)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_results()
