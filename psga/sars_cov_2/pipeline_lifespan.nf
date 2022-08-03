process pipeline_start {
  input:
    val metadata
    val run
    val sequencing_technology
    val kit

  output:

  script:
    session_id_file = "${PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH}/${run}_${workflow.sessionId}"

  """
  # save the command so that we can resume it using autoresumer.py, if needed
  echo "#!/bin/bash" > ${session_id_file}
  echo "export PSGA_OUTPUT_PATH=${PSGA_OUTPUT_PATH}" >> ${session_id_file}
  echo "nextflow run ${PSGA_ROOT_PATH}/psga -c ${PSGA_ROOT_PATH}/psga/sars_cov_2.config --run ${run} --sequencing_technology ${sequencing_technology} --metadata ${metadata} --kit ${kit} -resume \\\$1" >> ${session_id_file}
  """
}
