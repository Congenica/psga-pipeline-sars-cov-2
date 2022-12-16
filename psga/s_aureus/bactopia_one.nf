/*
 * Run: Run Bactopia analysis of s_aureus
 */
process bactopia_one {
  publishDir "${params.output_path}/bactopia_one", mode: 'copy', overwrite: true

  tag "${task.index} - ${ch_input_files}"

  input:
    path ch_input_files
    path ch_metadata_files

  output:
    path ch_input_files, emit: ch_input_files
    path "local-single-sample/bactopia/nf-reports/bactopia-trace.txt", emit: bactopia_trace
    path "local-single-sample/*/variants/*/*.tab", emit: variants_summary_tab
    path "local-single-sample/*/variants/*/*.bam", emit: bam_file
    path "local-single-sample/*/variants/*/*.bam.bai", emit: bam_bai_file
    path "local-single-sample/*/variants/*/*.vcf.gz", emit: all_variants_all_vcfs     // outputs all of
    path "local-single-sample/*/variants/*/*consensus.fa.gz", emit: consensus_fasta
    path "local-single-sample/*/variants/*/*consensus.subs.fa.gz", emit: consensus_fasta_subs
    path "local-single-sample/*/variants/*/*consensus.subs.masked.fa.gz", emit: consensus_fasta_masked
    path "local-single-sample/*/quality-control/summary/*_R1-final_fastqc.zip", emit: final_fastqc_r1
    path "local-single-sample/*/quality-control/summary/*_R2-final_fastqc.zip", emit: final_fastqc_r2
    path "local-single-sample/*/quality-control/summary/*_R1-original_fastqc.zip", emit: original_fastqc_r1
    path "local-single-sample/*/quality-control/summary/*_R2-original_fastqc.zip", emit: original_fastqc_r2
    path "local-single-sample/*/mlst/default/blast/*blast.json", emit: mlst_json
    path "local-single-sample/bactopia-tools/mlst/mlst/*/*.tsv", emit: mlst_tsv
    path "local-single-sample/bactopia-tools/mlst/mlst/*/logs/mlst/nf-mlst.err", emit: mlst_err
    path "local-single-sample/*/logs/assembly_qc/checkm.log", emit: checkm_log
    path "local-single-sample/*/assembly/checkm/bins/*/genes.gff", emit: checkm_genes_gff
    path "local-single-sample/*/assembly/checkm/bins/*/hmmer.tree.txt", emit: checkm_hmmer_tree
    path "local-single-sample/*/assembly/checkm/checkm-results.txt", emit: checkm_results_txt
    path "local-single-sample/logs/custom_dumpsoftwareversions/versions.yml", emit: software_versions
    path "local-single-sample/software_versions.yml", emit: all_software_versions
    path "local-single-sample/bactopia-tools/mykrobe/mykrobe/*/*.json", emit: mykrobe_json
    path "local-single-sample/bactopia-tools/mykrobe/mykrobe/*/logs/mykrobe/nf-mykrobe.log", emit: mykrobe_log
    path "local-single-sample/bactopia-tools/mykrobe/mykrobe/*/logs/mykrobe/versions.yml", emit: mykrobe_versions
    path "local-single-sample/*/antimicrobial-resistance/*-protein-report.txt", emit: antimicrobial_protein_report
    path "local-single-sample/*/antimicrobial-resistance/*-gene-report.txt", emit: antimicrobial_gene_report
    path "local-single-sample/*/annotation/*.tsv", emit: annotation_tsv
    path "local-single-sample/*/assembly/*.json", emit: assembly_json
    path "quast_assembly.zip", emit: quast_assembly_zip
    path "local-single-sample/*/variants/*/*.txt", emit: variants_txt_for_csv_file
    path "sample-*-csv-files/annotation-summary.txt", emit: annotation_summary
    path "local-single-sample/bactopia-tools/mlst/mlst/software_versions.yml", emit: mlst_versions
    path "genome_assembly_name.txt", emit: genome_assembly_name
    path "*-resultfiles.json", emit:ch_result_files_json
    path "*-results.csv", emit:ch_results_csv


  script:
  """
#!/usr/bin/env python
import os
import re
import shutil
import json

with open('hello.txt', 'w') as outfile:
    # Sample ID is not passed to script so get it from the sequence file names
    sample_id = None
    for line in os.listdir(os.getcwd()):
        match = re.match('(.+)_1.fastq.gz', line)
        if match:
            sample_id = match.group(1)
    if not sample_id:
        raise Exception('Cannot get sample_id from input filenames')

    bactopia_cmd = 'bactopia --R1 {sample_id}_1.fastq.gz --R2 {sample_id}_2.fastq.gz --sample {sample_id} --datasets /datasets/ --species "Staphylococcus aureus" --coverage 100 --genome_size median --outdir local-single-sample --max_cpus 10 -qs 1 --max_retry 1 --run_checkm'.format(
        sample_id=sample_id)

    outfile.write(bactopia_cmd + " ## ")
    outfile.write(str(os.getcwd()) + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

    # Run main bactopia command
    os.system(bactopia_cmd)

    # Run extra bactopia commands
    with open('includes.txt', 'w') as includes_file:
        includes_file.write(sample_id)
    os.system('bactopia --wf amrfinderplus --bactopia local-single-sample --organism Staphylococcus_aureus --include includes.txt')
    os.system('bactopia --wf mlst --bactopia local-single-sample --include includes.txt')
    os.system('bactopia --wf mykrobe --bactopia local-single-sample --mykrobe_species staph --include includes.txt')

    # Make an archive of the quast assembly so all the small files are in one place
    shutil.make_archive('quast_assembly', 'zip', 'local-single-sample/{sample_id}/assembly/quast'.format(sample_id=sample_id))

    # Copy the stuff for the results.csv prefixed with sample_id - avoids name clashes as 2 files are sample_id.txt
    # Directory name is odd because we need to match it without knowing sample_id
    outdir = 'sample-{sample_id}-csv-files'.format(sample_id=sample_id)
    os.mkdir(outdir)
    shutil.copy(
        'local-single-sample/{sample_id}/annotation/{sample_id}.txt'.format(sample_id=sample_id),
        '{outdir}/annotation-summary.txt'.format(outdir=outdir)
    )

    # TODO find a better way to get the genome assembly name.
    # This will find the genome assembly directory and copy it's name
    # Required because the resultfiles.json needs to know the directory names
    _ ,gcf_name, _ = next(os.walk('local-single-sample/{sample_id}/variants'.format(sample_id=sample_id), topdown=True))
    genome_assembly_name = gcf_name[0]
    print("genome_assembly_name: " + genome_assembly_name)
    with open('genome_assembly_name.txt', 'w') as genome_assembly_name_file:
        json.dump({"genome_assembly_name": genome_assembly_name, "sample_id": sample_id}, genome_assembly_name_file)

    outfile.write("RUN COMMAND" + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

### For some reason python can't find this script with a system call when run by nextflow. Just copy the script here
### Make this run from the standalone script

# flake8: noqa
from os.path import join as join_path  # used to join FS paths and S3 URIs
from pathlib import Path, PosixPath
from typing import Dict, List, Set
from functools import partial, reduce
import random
import csv
import json
from json import JSONEncoder
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



def generate_results(
    output_path: str,
) -> None:
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
    generate_results("${params.output_path}")



  """
}
