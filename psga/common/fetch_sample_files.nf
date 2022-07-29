/*
 * Stage the sample file
 */
process stage_sample_file {
  tag "${task.index} - ${file_1}"
  input:
    // use val for file_2 as it is part of the tuple, but is not used. This avoids the unwanted warning.
    tuple val(sample_id), path(file_1), val(file_2)

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
