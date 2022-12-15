process run_analysis {
  publishDir "${params.output_path}/analysis", mode: 'copy', overwrite: true, pattern: "*_analysis.txt"

  tag "${task.index} - ${ch_input_files}"

  input:
    path ch_input_files

  output:
    path ch_input_files, emit: ch_input_files
    path "*_analysis.txt", emit: ch_analysis_report

  shell:
  '''
  for fq in "!{ch_input_files}"; do
      sample_id="$(echo $fq | cut -d'.' -f1 | cut -d'_' -f1)"
      break
  done
  echo "Create dummy analysis report" > ${sample_id}_analysis.txt
  '''
}
