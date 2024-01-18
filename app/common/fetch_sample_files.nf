/*
 * Stage the sample file
 */
process stage_sample_file {
  tag "${task.index} - ${file_1}"
  input:
    tuple val(sample_id), path(file_1)

  output:
    path(file_1), emit: ch_sample_files

  script:
  """
  """
}

/*
 * Stage the sample file pair
 */
process stage_sample_file_pair {
  tag "${task.index} - [${file_1}, ${file_2}]"
  input:
    tuple val(sample_id), path(file_1), path(file_2)

  output:
    tuple path(file_1), path(file_2), emit: ch_sample_files

  script:
  """
  """
}
