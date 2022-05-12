/*
 * Run: fastqc
 */
process fastqc {
  tag "${task.index} - ${ch_input_files}"

  input:
    path ch_input_files

  output:
    path ch_input_files, emit: ch_input_files
    path "*_fastqc.html", emit: ch_fastqc_html_report
    path "*_fastqc.zip", emit: ch_fastqc_zip_report

  shell:
  '''
  # at this stage, all our sample files have extension: .fastq.gz (illumina) or .fastq (nanopore)
  for fq in `ls *.fastq*`; do
      fastqc ${fq}
  done
  '''
}

/*
 * Process to store all fastq_reports
 */
process store_fastqc_reports {
  publishDir "${PSGA_OUTPUT_PATH}/fastqc", mode: 'copy', overwrite: true

  input:
    path ch_fastqc_html_reports
    path ch_fastqc_zip_reports

  output:
    path ch_fastqc_html_reports, emit: ch_fastq_html_reports
    path ch_fastqc_zip_reports, emit: ch_fastq_zip_reports

  shell:
  '''
  '''
}
