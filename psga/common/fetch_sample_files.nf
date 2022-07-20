/*
 * Stage the sample file and check file integrity.
 * The returned file is renamed using the sample_id.
 * An error is raised if the integrity check fails.
 * Set the error strategy for this process if the pipeline should ignore this sample.
 */
process stage_sample_file {
  tag "${task.index} - ${staged_file_1}"
  input:
    val(analysis_run)
    val(file_extension)
    tuple val(sample_id), path(file_1)

  output:
    path(staged_file_1), emit: ch_sample_files

  script:
    staged_file_1 = "${sample_id}${file_extension}"
  """
  if [ ! -e ${staged_file_1} ] ; then
      # only create this symlink if the there isn't a file with name $staged_file_1.
      # This is only the case when the uploaded files are named after their lab names.
      # don't use ln -sf because nextflow removes both files.
      ln -s ${file_1} ${staged_file_1}
  fi
  """
}

/*
 * Stage the sample file pair and check file integrity.
 * The returned files are renamed using the sample_id.
 * An error is raised if the integrity check fails.
 * Set the error strategy for this process if the pipeline should ignore this sample.
 */
process stage_sample_file_pair {
    tag "${task.index} - [${staged_file_1}, ${staged_file_2}]"
    input:
      val(analysis_run)
      val(file_extension)
      tuple val(sample_id), path(file_1), path(file_2)

    output:
      tuple path(staged_file_1), path(staged_file_2), emit: ch_sample_files

    script:
      staged_file_1 = "${sample_id}_1${file_extension}"
      staged_file_2 = "${sample_id}_2${file_extension}"
    """
    if [ ! -e ${staged_file_1} ] ; then
        ln -s ${file_1} ${staged_file_1}
    fi
    if [ ! -e ${staged_file_2} ] ; then
        ln -s ${file_2} ${staged_file_2}
    fi
    """
}

workflow notifications {
    take:
        ch_files_matching_metadata
    main:

        ch_files_matching_metadata.ifEmpty {
            log.error """\
              ERROR: No file was found for samples in metadata.
                Aborting!
            """
            System.exit(1)
        }

    emit:
        ch_files_matching_metadata
}

workflow select_sample_file {
    take:
        ch_metadata
        file_extension
    main:

        ch_metadata
            .splitCsv(header:true)
            .map{ row-> tuple(row.sample_id, row.file_1) }
            .set { ch_sample_filepaths }
        // e.g. ch_sample_filepaths: [[sample_name1, fastq], ...]

        ch_input_files = Channel.empty()

        ch_sample_files = stage_sample_file(
            params.run,
            file_extension,
            ch_sample_filepaths
        )

        ch_sample_files
             .map { it -> [ it ] }
             .set{ ch_files_matching_metadata }

        notifications(ch_files_matching_metadata)

        ch_selected_sample_files = notifications.out.ch_files_matching_metadata

    emit:
        ch_selected_sample_files
}

workflow select_sample_file_pair {
    take:
        ch_metadata
        file_extension
    main:

        ch_metadata
            .splitCsv(header:true)
            .map{ row-> tuple(row.sample_id, row.file_1, row.file_2) }
            .set { ch_sample_filepaths }
        // e.g. ch_sample_filepaths: [[sample_name1, fastq_1, fastq_2], ...]

        ch_sample_files = stage_sample_file_pair(
            params.run,
            file_extension,
            ch_sample_filepaths
        )

        ch_sample_files
             .map { it -> [ it[0], it[1] ] }
             .set{ ch_files_matching_metadata }

        notifications(ch_files_matching_metadata)

        ch_selected_sample_files = notifications.out.ch_files_matching_metadata

    emit:
        ch_selected_sample_files
}
