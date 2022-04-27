include { makeBamSearchPath } from './utils.nf'
include { filter_samples_with_one_file } from './filtering.nf'
include { notifications } from './filtering.nf'


workflow filter_bam_matching_with_metadata{
    take:
        ch_all_samples_with_metadata_loaded
        ch_samples_with_metadata_loaded
        ch_qc_passed_samples
    main:
        ch_file_paths = makeBamSearchPath()

        Channel
            .fromPath(ch_file_paths)
            .map { file -> tuple(file.baseName, file) }
            .set{ ch_sample_files }

        filter_samples_with_one_file(
          ch_all_samples_with_metadata_loaded,
          ch_samples_with_metadata_loaded,
          ch_qc_passed_samples,
          ch_sample_files
        )

        notifications(
          filter_samples_with_one_file.out.ch_metadata_sample_mismatch_search,
          filter_samples_with_one_file.out.ch_mismatching_metadata_sample_names,
          filter_samples_with_one_file.out.ch_sample_qc_pass_matches_search,
          filter_samples_with_one_file.out.ch_files_matching_metadata
        )

        ch_selected_sample_files = notifications.out.ch_files_matching_metadata
        store_notification_missing_files = notifications.out.store_notification_missing_files

    emit:
        ch_selected_sample_files
        store_notification_missing_files
}
