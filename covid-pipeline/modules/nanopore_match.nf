include { makeNanoporeSearchPath } from './utils.nf'

include { append_metadata_match_to_sample_file }  from './utils.nf'
include { store_mismatching_files }  from './utils.nf'
include { append_match_to_values_list as append_match_to_samples }  from './utils.nf'
include { append_match_to_values_list as append_qc_pass_match_to_sample_files }  from './utils.nf'

include { store_notification_with_values_list as store_notification_missing_files } from './utils.nf'
include { store_notification_with_values_list as store_notification_updated } from './utils.nf'
include { store_notification_with_values_list as store_notification_sample_files_only } from './utils.nf'
include { store_notification_with_values_list as store_notification_processed_already } from './utils.nf'


workflow filter_nanopore_matching_with_metadata{
    take:
        ch_all_samples_with_metadata_loaded
        ch_samples_with_metadata_loaded
        ch_qc_passed_samples
        ch_updated_samples
    main:
        ch_file_paths = makeNanoporeSearchPath()

        Channel
            .fromPath(ch_file_paths)
            .map { file -> tuple(file.baseName, file) }
            .set{ ch_sample_files }

        ch_sample_files
            .map { it ->  it[0] }
            .flatten()
            .set{ ch_sample_names }

        // Matching our extracted sample names with samples loaded from metadata
        ch_sample_files_matches = append_metadata_match_to_sample_file(
            ch_sample_files,
            ch_all_samples_with_metadata_loaded.collect()
        )

        ch_sample_files_matches
            .branch { it ->
                matched: it[1]
                    return it
                mismatched: true
                    return it
            }
            .set{ ch_sample_files_grouped_by_metadata_search }

        ch_sample_files_grouped_by_metadata_search.matched
             .map { it -> it[1] }
             .flatten()
             .set{ ch_files_matching_metadata }

        ch_sample_files_grouped_by_metadata_search.mismatched
             .map { it -> it[1] }
             .flatten()
             .set{ ch_files_mismatching_metadata }
        ch_sample_files_grouped_by_metadata_search.mismatched
             .map { it -> it[0] }
             .flatten()
             .set{ ch_mismatching_metadata_sample_names }

        store_mismatching_files(ch_files_mismatching_metadata)

        ch_metadata_sample_matches = append_match_to_samples(
            ch_samples_with_metadata_loaded,
            ch_sample_names.collect()
        )

        ch_metadata_sample_matches
            .branch { it ->
                matched: it[1]
                    return it[0]
                mismatched: true
                    return it[0]
            }
            .set{ ch_metadata_sample_matches_search }

        ch_sample_qc_pass_matches = append_qc_pass_match_to_sample_files(
            ch_sample_names,
            ch_qc_passed_samples.collect()
        )

        ch_sample_qc_pass_matches
            .branch { it ->
                matched: it[1]
                    return it[0]
                mismatched: true
                    return it[0]
            }
            .set{ ch_sample_qc_pass_matches_search }

        store_notification_missing_files(
          "${workflow.start}-samples_missing_files.txt",
          ch_metadata_sample_matches_search.mismatched.collect()
        ).subscribe onNext: {
            log.error("ERROR: Found samples with missing files. See notification file: samples_missing_files.txt. Abort!")
            System.exit(1)
        }

        store_notification_updated(
          "${workflow.start}-updated_samples.txt",
          ch_updated_samples.collect()
        ).subscribe onNext: {
            log.warn("Found samples with updated metadata. See notification file: updated_samples.txt.")
        }

        store_notification_sample_files_only(
          "${workflow.start}-samples_missing_metadata.txt",
          ch_mismatching_metadata_sample_names.collect()
        ).subscribe onNext: {
            log.warn("Found sample files missing metadata. These will be skipped. See notification file: samples_missing_metadata.txt.")
        }

        store_notification_processed_already(
          "${workflow.start}-samples_already_processed.txt",
          ch_sample_qc_pass_matches_search.matched.collect()
        ).subscribe onNext: {
            log.warn("Found samples already marked as QC_PASS by ncov2019-artic. They will be re-processed and results will be overwritten. See notification file: samples_already_processed.txt.")
        }

        ch_files_matching_metadata.ifEmpty {
            log.error """\
              ERROR: No file was found for samples in metadata.
                This may be caused by a failure in loading samples from the metadata file to the database.
                Aborting!
            """
            System.exit(1)
        }

    emit:
        ch_files_matching_metadata
        store_notification_missing_files
}
