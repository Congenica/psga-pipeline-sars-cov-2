/*
 * Check the metadata in the database
 */
process check_metadata {
  publishDir "${PSGA_OUTPUT_PATH}/notifications", mode: 'copy', overwrite: true, pattern: 'samples_with_{invalid,valid}_metadata.txt'
  publishDir "${PSGA_OUTPUT_PATH}/logs", mode: 'copy', overwrite: true, pattern: '*.log'

  tag "${metadata}"

  input:
    path metadata
    val analysis_run_name
    val filetype
    val ncov_workflow

  output:
    path metadata, emit: ch_metadata
    path "check_metadata.done", emit: ch_metadata_checked
    path "samples_with_invalid_metadata.txt", emit: ch_samples_with_invalid_metadata_file
    path "samples_with_valid_metadata.txt", emit: ch_samples_with_valid_metadata_file
    path "*.log"

  shell:
  '''
  # convert nextflow variables to Bash so that the same format is used

  metadata=!{metadata}
  analysis_run_name=!{analysis_run_name}
  filetype=!{filetype}
  ncov_workflow=!{ncov_workflow}

  metadata_checked="check_metadata.done"
  samples_with_invalid_metadata_file="samples_with_invalid_metadata.txt"
  samples_with_valid_metadata_file="samples_with_valid_metadata.txt"

  python ${PSGA_ROOT_PATH}/scripts/sars_cov_2/check_metadata.py \
    --metadata-path "${metadata}" \
    --analysis-run-name "${analysis_run_name}" \
    --input-file-type "${filetype}" \
    --ncov-workflow "${ncov_workflow}" \
    --samples-with-invalid-metadata-file "${samples_with_invalid_metadata_file}" \
    --samples-with-valid-metadata-file "${samples_with_valid_metadata_file}"

  touch ${metadata_checked}
  '''
}