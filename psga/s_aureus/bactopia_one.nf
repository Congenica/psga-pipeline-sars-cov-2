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
    path "hello.txt", emit: hello
    path "metadata.csv", emit: md_out
    path "local-single-sample/bactopia/nf-reports/bactopia-trace.txt", emit: bactopia_trace
//     path "local-single-sample/*/variants/*/*.bam", emit: bam_file
//     path "local-single-sample/*/variants/*/*annotated.vcf.gz", emit: annotated_vcf
//     path "local-single-sample/*/variants/*/*consensus.fa.gz", emit: consensus_fasta
//     path "local-single-sample/9729bce7-f0a9-4617-b6e0-6145307741d1/logs/call_variants/refseq-genomes/versions.yml", emit: call_variants_versions
//     path "local-single-sample/9729bce7-f0a9-4617-b6e0-6145307741d1/logs/antimicrobial-resistance/9729bce7-f0a9-4617-b6e0-6145307741d1-protein-report.txt", emit: antimicrobial_resistance
    // TODO work out how to get sample name in this part of the file ... or call bactopia with a constant sample name
    stdout

  script:
  """
#!/usr/bin/env python
import csv
import os
import re

with open('hello.txt', 'w') as outfile:
    # Sample ID is not passed to script so get it from the sequence file names
    sample_id = None
    for line in os.listdir(os.getcwd()):
        match = re.match('(.+)_1.fastq.gz', line)
        if match:
            sample_id = match.group(1)
    if not sample_id:
        raise Exception('Cannot get sample_id from input filenames')

    bactopia_cmd = 'bactopia --R1 {sample_id}_1.fastq.gz --R2 {sample_id}_2.fastq.gz --sample {sample_id} --datasets /datasets/ --species "Staphylococcus aureus" --coverage 100 --genome_size median --outdir local-single-sample --max_cpus 10 -qs 1 --max_retry 1'.format(
        sample_id=sample_id)

    outfile.write(bactopia_cmd + " ## ")
    outfile.write(str(os.getcwd()) + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

    os.system(bactopia_cmd)

    outfile.write("RUN COMMAND" + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

  """
}
