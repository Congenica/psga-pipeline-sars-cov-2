process pipeline_started {
  input:
    val run
    val pipeline_workflow
    val filetype
    val scheme_repo_url
    val scheme_dir
    val scheme
    val scheme_version

  output:

  script:
    session_id_file = "${PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH}/${run}_${workflow.sessionId}"

  """
  # save the command so that we can resume it using autoresumer.py, if needed
  echo "#!/bin/bash" > ${session_id_file}
  echo "export PSGA_INPUT_PATH=${PSGA_INPUT_PATH}" >> ${session_id_file}
  echo "export PSGA_OUTPUT_PATH=${PSGA_OUTPUT_PATH}" >> ${session_id_file}
  echo "nextflow run ${PSGA_ROOT_PATH}/psga --run ${run} --workflow ${pipeline_workflow} --filetype ${filetype} --scheme_repo_url ${scheme_repo_url} --scheme_dir ${scheme_dir} --scheme ${scheme} --scheme_version ${scheme_version} -resume \\\$1" >> ${session_id_file}
  """
}

process pipeline_completed {
  input:
    val run
    file ncov_qc_sample_submitted_complete
    file pangolin_sample_submitted_complete

  output:

  shell:
  '''
  # use shell because we need a wildcard
  session_id_file="!{PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH}/!{run}_!{workflow.sessionId}"
  # remove file and any attempts as no longer needed
  rm ${session_id_file}*
  '''
}