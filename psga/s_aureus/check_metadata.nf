/*
 * Check the metadata in the database
 */
process check_metadata {

  tag "${metadata}"

  input:
    path metadata

  output:
    path metadata, emit: ch_metadata

  shell:
  '''
  metadata=!{metadata}

  metadata_checked="check_metadata.done"
  touch ${metadata_checked}
  '''
}
