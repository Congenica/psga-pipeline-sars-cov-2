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
      tuple val(sample_id), path(file_1), val(md5_1)

    output:
      path(staged_file_1), emit: ch_sample_files

    script:
      staged_file_1 = "${sample_id}${file_extension}"
    """
    ln -s ${file_1} ${staged_file_1}

    python ${PSGA_ROOT_PATH}/scripts/validation/check_file_integrity.py --analysis-run-name ${analysis_run} --sample-name ${sample_id} --input-path ${staged_file_1} --expected-md5 ${md5_1}
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
      tuple val(sample_id), path(file_1), path(file_2), val(md5_1), val(md5_2)

    output:
      tuple path(staged_file_1), path(staged_file_2), emit: ch_sample_files

    script:
      staged_file_1 = "${sample_id}_1${file_extension}"
      staged_file_2 = "${sample_id}_2${file_extension}"
    """
    ln -s ${file_1} ${staged_file_1}
    ln -s ${file_2} ${staged_file_2}

    python ${PSGA_ROOT_PATH}/scripts/validation/check_file_integrity.py --analysis-run-name ${analysis_run} --sample-name ${sample_id} --input-path ${staged_file_1} --expected-md5 ${md5_1}
    python ${PSGA_ROOT_PATH}/scripts/validation/check_file_integrity.py --analysis-run-name ${analysis_run} --sample-name ${sample_id} --input-path ${staged_file_2} --expected-md5 ${md5_2}
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
            .map{ row-> tuple(row.sample_id, row.file_1, row.md5_1) }
            .set { ch_sample_filepaths }
        // e.g. ch_sample_filepaths: [[sample_name1, fastq, md5], ...]

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
            .map{ row-> tuple(row.sample_id, row.file_1, row.file_2, row.md5_1, row.md5_2) }
            .set { ch_sample_filepaths }
        // e.g. ch_sample_filepaths: [[sample_name1, fastq_1, fastq_2, md5_1, md5_2], ...]

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
