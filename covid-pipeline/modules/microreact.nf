/*
 * Prepare input tsv for microreact
 */
process prepare_microreact_tsv {
  publishDir "${COVID_PIPELINE_MICROREACT_PATH}/${workflow.sessionId}", mode: 'copy', overwrite: true
  publishDir "${COVID_PIPELINE_MICROREACT_PATH}/latest", mode: 'copy', overwrite: true

  input:
    path ncov_qc_to_db_submit_completion_flag
    path pangolin_to_db_submit_completion_flag

  output:
    path microreact_file

  script:
    microreact_file = params.microreact_tsv

  """
  python /app/scripts/generate_microreact_input.py --output ${microreact_file}
  """
}
