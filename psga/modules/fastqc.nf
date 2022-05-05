/*
 * Run: fastqc
 */
process fastqc {
  tag "${task.index} - ${fastq_file}"

  input:
    path fastq_file

  output:
    path "fastqc.done", emit: ch_fastqc_done
    path "*_fastqc.html", emit: ch_fastqc_html_report
    path "*_fastqc.zip", emit: ch_fastqc_zip_report

  shell:
  '''
  # at this stage, all our sample files have extension: .fastq.gz (illumina) or .fastq (nanopore)
  for fq in `ls *.fastq*`; do
      fastqc ${fq}
  done

  touch "fastqc.done"
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
