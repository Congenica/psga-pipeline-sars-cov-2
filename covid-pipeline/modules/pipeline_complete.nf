/*
 * Delete the pipeline_started file to mark this pipeline run as complete.
 * This step will not return an error if removing the file fails, for example if it didn't exist.
 * This is to prevent errors if the pipeline is run manually.
 */
process pipeline_complete {
  input:
    file ncov_qc_sample_submitted_complete
    file pangolin_sample_submitted_complete

  output:

  script:
    pipeline_started_file = "${COVID_PIPELINE_INPUT_PATH}/" + "pipeline_started"
    pipeline_complete_file = "${COVID_PIPELINE_INPUT_PATH}/" + "pipeline_complete"

  """
  touch ${pipeline_complete_file}
  rm ${pipeline_started_file} || true
  """
}
