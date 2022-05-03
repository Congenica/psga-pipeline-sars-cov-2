/*
 * Check the metadata in the database
 */
process check_metadata {
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
    path "current_samples_with_metadata.txt", emit: ch_current_session_samples_with_metadata_file

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
  current_samples_with_metadata="current_samples_with_metadata.txt"

  touch ${current_samples_with_metadata}

  [[ "${load_missing_samples}" == "true" ]] && load_samples_flag="--load-missing-samples" || load_samples_flag=""

  python ${PSGA_ROOT_PATH}/scripts/check_metadata.py \
    --file "${metadata}" \
    --analysis-run-name "${analysis_run_name}" \
    --primer-scheme-name "${scheme}" \
    --primer-scheme-version "${scheme_version}" \
    --input-file-type "${filetype}" \
    --workflow "${ncov_workflow}" \
    --output-current-samples-with-metadata "${current_samples_with_metadata}" \
    --pipeline-version "${pipeline_version}" \
    ${load_samples_flag}

  touch ${metadata_checked}
  '''
}
