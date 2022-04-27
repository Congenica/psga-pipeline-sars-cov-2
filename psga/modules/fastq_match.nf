include { makeFastqSearchPath } from './utils.nf'
include { filter_samples_with_two_files } from './filtering.nf'
include { notifications } from './filtering.nf'


workflow filter_fastq_matching_with_metadata{
    take:
        ch_all_samples_with_metadata_loaded
        ch_samples_with_metadata_loaded
        ch_qc_passed_samples
    main:
        // Original usage of grouping illumina fastq file pairs and extracting sample name
        ch_file_paths = makeFastqSearchPath(
            params.illuminaPrefixes,
            params.illuminaSuffixes,
            params.fastq_exts
        )
        Channel.fromFilePairs( ch_file_paths, flat: true)
            .filter{ !( it[0] =~ /Undetermined/ ) }
            .set{ ch_sample_files }
        // e.g. ch_sample_files: [[sample_name1, fastq1_1, fastq1_2], ...]

        filter_samples_with_two_files(
          ch_all_samples_with_metadata_loaded,
          ch_samples_with_metadata_loaded,
          ch_qc_passed_samples,
          ch_sample_files
        )

        notifications(
          filter_samples_with_two_files.out.ch_metadata_sample_mismatch_search,
          filter_samples_with_two_files.out.ch_mismatching_metadata_sample_names,
          filter_samples_with_two_files.out.ch_sample_qc_pass_matches_search,
          filter_samples_with_two_files.out.ch_files_matching_metadata
        )

        ch_selected_sample_files = notifications.out.ch_files_matching_metadata
        store_notification_missing_files = notifications.out.store_notification_missing_files

    emit:
        ch_selected_sample_files
        store_notification_missing_files
}
