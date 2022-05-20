process pipeline_end {
  input:
    val run
    path ch_sample_results_submitted

  output:

  shell:
  '''
  # use shell because we need a wildcard
  session_id_file="!{PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH}/!{run}_!{workflow.sessionId}"
  # remove file and any attempts as no longer needed
  rm -f ${session_id_file}*
  '''
}
