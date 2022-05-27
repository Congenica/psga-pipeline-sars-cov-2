/*
 * Run: fastqc
 */
process fastqc {
  publishDir "${PSGA_OUTPUT_PATH}/fastqc", mode: 'copy', overwrite: true, pattern: '*_fastqc.{html,zip}'

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
