/*
 * Load the metadata into the database
 */
process load_metadata {
  input:
    path ch_metadata_tsv_file
    val ch_analysis_run_name

  output:
    path metadata_load_done, emit: ch_metadata_load_done
    path all_samples_with_metadata, emit: ch_all_samples_with_metadata_file
    path current_samples_with_metadata, emit: ch_current_session_samples_with_metadata_file
    path samples_with_qc_pass, emit: ch_all_samples_ncov2019_artic_qc_passed_file
    path updated_samples, emit: ch_current_session_updated_samples_file

  script:
    metadata_load_done = "load_metadata.done"
    all_samples_with_metadata = "all_samples_with_metadata.txt"
    current_samples_with_metadata = "current_samples_with_metadata.txt"
    samples_with_qc_pass = "samples_with_qc_pass.txt"
    updated_samples = "updated_samples.txt"

  """
  touch ${all_samples_with_metadata}
  touch ${current_samples_with_metadata}
  touch ${samples_with_qc_pass}
  touch ${updated_samples}

  python /app/scripts/load_metadata_to_db.py \
    --file "${ch_metadata_tsv_file}" \
    --analysis-run-name "${ch_analysis_run_name}" \
    --output-all-samples-with-metadata "${all_samples_with_metadata}" \
    --output-current-samples-with-metadata "${current_samples_with_metadata}" \
    --output-samples-with-qc-pass "${samples_with_qc_pass}" \
    --output-samples-updated "${updated_samples}" \
    --pipeline-version "${workflow.manifest.version}"
  touch ${metadata_load_done}
  """
}
