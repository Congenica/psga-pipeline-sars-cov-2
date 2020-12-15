/*
 * Load the I-SEHA metadata into the database
 */
process load_iseha_metadata {
  input:
    path ch_iseha_metadata_tsv_file

  output:
    path iseha_metadata_load_done, emit: ch_iseha_metadata_load_done
    path samples_with_metadata, emit: ch_samples_with_metadata_file

  script:
    iseha_metadata_load_done = "load_iseha_metadata.done"
    samples_with_metadata = "all_samples_with_metadata.txt"

  """
  python /app/scripts/load_iseha_metadata.py \
    --file "${ch_iseha_metadata_tsv_file}" \
    --output_file_for_samples_with_metadata "${samples_with_metadata}"
  touch ${iseha_metadata_load_done}
  """
}
