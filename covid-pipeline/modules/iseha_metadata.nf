/*
 * Load the I-SEHA metadata into the database
 */
process load_iseha_metadata {
  input:
    path ch_iseha_metadata_tsv_file

  output:
    path iseha_metadata_load_done, emit: ch_iseha_metadata_load_done
    path all_samples_with_metadata, emit: ch_all_samples_with_metadata_file
    path current_samples_with_metadata, emit: ch_current_session_samples_with_metadata_file
    path samples_with_qc_pass, emit: ch_all_samples_ncov2019_artic_qc_passed_file
    path updated_samples, emit: ch_current_session_updated_samples_file

  script:
    iseha_metadata_load_done = "load_iseha_metadata.done"
    all_samples_with_metadata = "all_samples_with_metadata.txt"
    current_samples_with_metadata = "current_samples_with_metadata.txt"
    samples_with_qc_pass = "samples_with_qc_pass.txt"
    updated_samples = "updated_samples.txt"

  """
  touch ${all_samples_with_metadata}
  touch ${current_samples_with_metadata}
  touch ${samples_with_qc_pass}
  touch ${updated_samples}

  python /app/scripts/load_iseha_metadata.py \
    --file "${ch_iseha_metadata_tsv_file}" \
    --output_all_samples_with_metadata "${all_samples_with_metadata}" \
    --output_current_samples_with_metadata "${current_samples_with_metadata}" \
    --output_samples_with_qc_pass "${samples_with_qc_pass}" \
    --output_samples_updated "${updated_samples}"
  touch ${iseha_metadata_load_done}
  """
}
