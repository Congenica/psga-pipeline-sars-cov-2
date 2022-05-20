/*
 * Check the metadata in the database
 */
process check_metadata {
  tag "${metadata}"

  input:
    val load_missing_samples
    path metadata
    val analysis_run_name
    val scheme
    val scheme_version
    val filetype
    val ncov_workflow

  output:
    path metadata, emit: ch_metadata
    path "check_metadata.done", emit: ch_metadata_checked
    path "samples_with_invalid_metadata.txt", emit: ch_samples_with_invalid_metadata_file
    path "samples_with_valid_metadata.txt", emit: ch_samples_with_valid_metadata_file

  shell:
  '''
  # convert nextflow variables to Bash so that the same format is used

  load_missing_samples=!{load_missing_samples}
  metadata=!{metadata}
  analysis_run_name=!{analysis_run_name}
  scheme=!{scheme}
  scheme_version=!{scheme_version}
  filetype=!{filetype}
  ncov_workflow=!{ncov_workflow}
  pipeline_version=!{workflow.manifest.version}

  metadata_checked="check_metadata.done"
  samples_with_invalid_metadata_file="samples_with_invalid_metadata.txt"
  samples_with_valid_metadata_file="samples_with_valid_metadata.txt"

  [[ "${load_missing_samples}" == "true" ]] && load_samples_flag="--load-missing-samples" || load_samples_flag=""

  python ${PSGA_ROOT_PATH}/scripts/sars_cov_2/check_metadata.py \
    --metadata-path "${metadata}" \
    --analysis-run-name "${analysis_run_name}" \
    --primer-scheme-name "${scheme}" \
    --primer-scheme-version "${scheme_version}" \
    --input-file-type "${filetype}" \
    --ncov-workflow "${ncov_workflow}" \
    --samples-with-invalid-metadata-file "${samples_with_invalid_metadata_file}" \
    --samples-with-valid-metadata-file "${samples_with_valid_metadata_file}" \
    --pipeline-version "${pipeline_version}" \
    ${load_samples_flag}

  touch ${metadata_checked}
  '''
}
