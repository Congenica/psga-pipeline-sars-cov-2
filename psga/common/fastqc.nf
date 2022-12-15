/*
 * Run: fastqc
 */
process fastqc {
  publishDir "${params.output_path}/fastqc", mode: 'copy', overwrite: true, pattern: '*_fastqc.zip'

  tag "${task.index} - ${ch_input_files}"

  input:
    path ch_input_files

  output:
    path ch_input_files, emit: ch_input_files
    path "*_fastqc.zip", emit: ch_fastqc_zip_report

  shell:
  '''
  for fq in "!{ch_input_files}"; do
      fastqc --limits !{params.fastqc_limits} ${fq}
  done
  '''
}
