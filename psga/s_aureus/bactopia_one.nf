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
//     path "hello.txt", emit: hello
    path "metadata.csv", emit: md_out
    path "local-single-sample/bactopia/nf-reports/bactopia-trace.txt", emit: bactopia_trace
    path "local-single-sample/*/variants/*/*.bam", emit: bam_file
    path "local-single-sample/*/variants/*/*.vcf.gz", emit: annotated_vcf
    path "local-single-sample/*/variants/*/*consensus.fa.gz", emit: consensus_fasta
    path "local-single-sample/*/quality-control/summary/*final_fastqc.zip", emit: final_fastqc
    path "local-single-sample/*/mlst/default/blast/*blast.json", emit: mlst_json
    path "local-single-sample/*/logs/assembly_qc/checkm.log", emit: checkm_log
    path "local-single-sample/*/assembly/checkm/bins/*/genes.gff", emit: checkm_genes_gff
    path "local-single-sample/*/assembly/checkm/bins/*/hmmer.tree.txt", emit: checkm_hmmer_tree
    path "local-single-sample/*/assembly/checkm/checkm-results.txt", emit: checkm_results_txt
    path "local-single-sample/*/antimicrobial-resistance/*-report.txt", emit: antimicrobial_resistance
    path "local-single-sample/logs/custom_dumpsoftwareversions/versions.yml", emit: software_versions

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

    bactopia_cmd = 'bactopia --R1 {sample_id}_1.fastq.gz --R2 {sample_id}_2.fastq.gz --sample {sample_id} --datasets /datasets/ --species "Staphylococcus aureus" --coverage 100 --genome_size median --outdir local-single-sample --max_cpus 10 -qs 1 --max_retry 1 --run_checkm'.format(
        sample_id=sample_id)

    outfile.write(bactopia_cmd + " ## ")
    outfile.write(str(os.getcwd()) + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

    os.system(bactopia_cmd)

    outfile.write("RUN COMMAND" + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

  """
}
