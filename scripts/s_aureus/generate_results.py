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
    #      TODO all local-single-sample/*/variants/*/*.vcf.gz -- problem with filename collisions
]


def build_csv_list(sample):
    """
    Gather all info for the CSV file
    Make a defined list of files that we are expecting.
    """
    CsvItem = namedtuple('CsvItem', 'header value')
    csv_list = list()
    # Add items to the list
    csv_list.append(CsvItem('Sample_Id', sample['SAMPLE_ID']))
    csv_list.append(CsvItem('QC_STATUS', 'PASSED'))  # For the demo everything passes

    with open('software_versions.yml') as software_versions:
        versions = list(csv.reader(software_versions, delimiter=":"))
        for row in versions:
            if 'checkm' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))
            if 'bactopia' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))
            if 'fastqc' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))

    with open('mlst_versions.yml') as mlst_versions:
        m_versions = list(csv.reader(mlst_versions, delimiter=":"))
        for row in m_versions:
            if 'mlst' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))

    with open('mykrobe_versions.yml') as mykrobe_versions:
        myk_versions = list(csv.reader(mykrobe_versions, delimiter=":"))
        for row in myk_versions:
            if 'mykrobe' in row[0]:
                csv_list.append(CsvItem('{}-version'.format(row[0].rstrip().lstrip()), row[1].lstrip().rstrip()))

    with open('annotation-summary.txt') as annotation_summary:
        annotations = csv.reader(annotation_summary, delimiter =":")
        for row in annotations:
            csv_list.append(CsvItem('Annotation-{}'.format(row[0]), row[1].lstrip().rstrip()))

    # Variant summary
    with open('{}.txt'.format(sample['SAMPLE_ID'])) as variant_summary_file:
        variant_summary = csv.reader(variant_summary_file, delimiter="\t")
        for row in variant_summary:
            if 'Variant' in row[0]:
                csv_list.append(CsvItem(row[0].rstrip().lstrip(), row[1].lstrip().rstrip()))

    with open('checkm-results.txt') as checkm_results_file:
        checkm_results = next(csv.DictReader(checkm_results_file, delimiter="\t"))
        csv_list.append(CsvItem('Marker lineage', checkm_results['Marker lineage']))
        csv_list.append(CsvItem('Completeness', checkm_results['Completeness']))
        csv_list.append(CsvItem('Contamination', checkm_results['Contamination']))
        csv_list.append(CsvItem('Strain heterogeneity', checkm_results['Strain heterogeneity']))

    with open('assembly.json') as assembly_results_file:
        assembly_results = json.loads(assembly_results_file.read())
        for key, value in assembly_results.items():
            if key == 'fasta':
                continue
            csv_list.append(CsvItem(key , value))

    with open('{sample_id}-protein-report.txt'.format(sample_id=sample['SAMPLE_ID'])) as antimicrobial_protein_file:
        resistance_genes = list()
        for row in csv.DictReader(antimicrobial_protein_file, delimiter="\t"):
            resistance_genes.append(row['Gene symbol'])
        csv_list.append(CsvItem('Antimicrobial resistance genes', ' '.join(resistance_genes)))

    return csv_list



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
    # TODO change to single sample -- get sample another way
    metadata = list(csv.DictReader(open(Path(metadata_file))))

    sample = metadata[0]

    csv_list = build_csv_list(sample)

    with open(output_csv_file, 'w', encoding='UTF8') as csv_out:
        writer = csv.writer(csv_out)
        writer.writerow([x.header for x in csv_list])
        writer.writerow([x.value for x in csv_list])

    # TODDO make this explicitly single sample. I will collate the results in another step
    # Build resultfiles.json
    json_out = dict()
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
