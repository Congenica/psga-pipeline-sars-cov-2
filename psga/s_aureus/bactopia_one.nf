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
    path "local-single-sample/*/variants/GCF_004153365/*.tab", emit: variants_summary_tab
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
    path "local-single-sample/*/variants/GCF_004153365/*.txt", emit: variants_txt_for_csv_file
    path "sample-*-csv-files/annotation-summary.txt", emit: annotation_summary

  script:
  """
#!/usr/bin/env python
import os
import re
import shutil

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

    outfile.write("RUN COMMAND" + " ## ")
    outfile.write(str(os.listdir(os.getcwd())) + " ## ")

  """
}
