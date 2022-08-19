process pipeline_start {
  input:

  output:

  script:
    session_id_file = "${PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH}/${params.run}_${workflow.sessionId}"

  shell:
  """
  # save the command so that we can resume it using autoresumer.py, if needed
  echo "#!/bin/bash" > ${session_id_file}
  echo "nextflow run ${PSGA_ROOT_PATH}/psga --run ${params.run} --sequencing_technology ${params.sequencing_technology} --metadata ${params.metadata} --kit ${params.kit} --output_path ${params.output_path} -resume \\\$1" >> ${session_id_file}
  """
}
