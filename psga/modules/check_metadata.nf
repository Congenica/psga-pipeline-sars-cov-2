/*
 * Check the metadata in the database
 */
process check_metadata {
  input:
    val load_missing_samples
    path metadata_file
    val analysis_run_name
    val scheme
    val scheme_version
    val filetype
    val ncov_workflow

  output:
    path "check_metadata.done", emit: ch_metadata_checked
    path "all_samples_with_metadata.txt", emit: ch_all_samples_with_metadata_file
    path "current_samples_with_metadata.txt", emit: ch_current_session_samples_with_metadata_file
    path "samples_with_qc_pass.txt", emit: ch_all_samples_ncov2019_artic_qc_passed_file

  shell:
  '''
  # convert nextflow variables to Bash so that the same format is used

  load_missing_samples=!{load_missing_samples}
  metadata_file=!{metadata_file}
  analysis_run_name=!{analysis_run_name}
  scheme=!{scheme}
  scheme_version=!{scheme_version}
  filetype=!{filetype}
  ncov_workflow=!{ncov_workflow}
  pipeline_version=!{workflow.manifest.version}

  metadata_checked="check_metadata.done"
  all_samples_with_metadata="all_samples_with_metadata.txt"
  current_samples_with_metadata="current_samples_with_metadata.txt"
  samples_with_qc_pass="samples_with_qc_pass.txt"

  touch ${all_samples_with_metadata}
  touch ${current_samples_with_metadata}
  touch ${samples_with_qc_pass}

  [[ "${load_missing_samples}" == "true" ]] && load_samples_flag="--load-missing-samples" || load_samples_flag=""

  python /app/scripts/check_metadata.py \
    --file "${metadata_file}" \
    --analysis-run-name "${analysis_run_name}" \
    --primer-scheme-name "${scheme}" \
    --primer-scheme-version "${scheme_version}" \
    --input-file-type "${filetype}" \
    --workflow "${ncov_workflow}" \
    --output-all-samples-with-metadata "${all_samples_with_metadata}" \
    --output-current-samples-with-metadata "${current_samples_with_metadata}" \
    --output-samples-with-qc-pass "${samples_with_qc_pass}" \
    --pipeline-version "${pipeline_version}" \
    ${load_samples_flag}

  touch ${metadata_checked}
  '''
}
