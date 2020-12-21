/*
 * Delete the pipeline_started file to mark this pipeline run as complete.
 * This step will not return an error if removing the file fails, for example if it didn't exist.
 * This is to prevent errors if the pipeline is run manually.
 */
process pipeline_complete {
  input:
    file strain_level_and_global_context_complete
    file strain_first_seen_complete
    file strain_prevalence_complete
    file sample_dump_complete

  output:

  script:
    pipeline_started_file = "${COVID_PIPELINE_FASTQ_PATH}/" + "pipeline_started"
    pipeline_complete_file = "${COVID_PIPELINE_FASTQ_PATH}/" + "pipeline_complete"

  """
  touch ${pipeline_complete_file}
  rm ${pipeline_started_file} || true
  """
}