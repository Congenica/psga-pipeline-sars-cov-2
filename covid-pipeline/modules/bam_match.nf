// Imports from 3rd party pipelines
include { makeBamSearchPath } from './utils.nf'

include { append_match_to_values_list as append_match_to_current_session_samples }  from './utils.nf'
include { append_match_to_values_list as append_qc_pass_match_to_bam_samples }  from './utils.nf'

include { store_notification_with_values_list as store_notification_no_bam_file } from './utils.nf'
include { store_notification_with_values_list as store_notification_updated } from './utils.nf'
include { store_notification_with_values_list as store_notification_bam_only } from './utils.nf'
include { store_notification_with_values_list as store_notification_processed_already } from './utils.nf'


workflow filter_bam_matching_with_metadata{
    take:
        ch_all_samples_with_metadata_loaded
        ch_current_session_samples_with_metadata_loaded
        ch_qc_passed_samples
        ch_updated_samples
    main:
        ch_bam_search_path = makeBamSearchPath()

        Channel
            .fromPath(ch_bam_search_path)
            .map { file -> tuple(file.baseName, file) }
            .set{ ch_bam_files }

        ch_bam_files
            .map { it ->  it[0] }
            .flatten()
            .set{ ch_bam_sample_names }

        // Matching our extracted sample names with samples loaded from I-SEHA
        ch_bam_file_pair_matches = append_metadata_match_to_sample_file(
            ch_bam_files,
            ch_all_samples_with_metadata_loaded.collect()
        )

        ch_bam_file_pair_matches
            .branch { it ->
                matched: it[1]
                    return it
                mismatched: true
                    return it
            }
            .set{ ch_bam_file_pairs_grouped_by_metadata_search }

        ch_bam_file_pairs_grouped_by_metadata_search.matched
             .map { it -> it[1] }
             .flatten()
             .set{ ch_bam_matching_metadata }

        ch_bam_file_pairs_grouped_by_metadata_search.mismatched
             .map { it -> it[1] }
             .flatten()
             .set{ ch_bam_mismatching_metadata }
        ch_bam_file_pairs_grouped_by_metadata_search.mismatched
             .map { it -> it[0] }
             .flatten()
             .set{ ch_mismatching_metadata_sample_names }

        store_bams_not_matched(ch_bam_mismatching_metadata)

        ch_current_session_metadata_sample_matches = append_match_to_current_session_samples(
            ch_current_session_samples_with_metadata_loaded,
            ch_bam_sample_names.collect()
        )

        ch_current_session_metadata_sample_matches
            .branch { it ->
                matched: it[1]
                    return it[0]
                mismatched: true
                    return it[0]
            }
            .set{ ch_current_session_metadata_sample_matches_search }

        ch_bam_sample_qc_pass_matches = append_qc_pass_match_to_bam_samples(
            ch_bam_sample_names,
            ch_qc_passed_samples.collect()
        )

        ch_bam_sample_qc_pass_matches
            .branch { it ->
                matched: it[1]
                    return it[0]
                mismatched: true
                    return it[0]
            }
            .set{ ch_bam_sample_qc_pass_matches_search }

        ch_current_session_metadata_sample_matches_search.mismatched.map {
            it -> log.warn """Sample ${it}, read from metadata file, has no matching bam file in current session. Sample will not be processed further"""
        }
        store_notification_no_bam_file(
          "${workflow.start}-no_bam_found_for_provided_sample.txt",
          ch_current_session_metadata_sample_matches_search.mismatched.collect()
        )

        ch_updated_samples.map {
            it -> log.info """Sample ${it} was updated with provided metadata"""
        }
        store_notification_updated(
          "${workflow.start}-updated_samples.txt",
          ch_updated_samples.collect()
        )

        ch_bam_file_pairs_grouped_by_metadata_search.mismatched.map {
            it -> log.warn """Sample ${it[0]} was not found in database with metadata.
                The following file will not be processed: ${it[1]}
            """
        }
        store_notification_bam_only(
          "${workflow.start}-no_metadata_provided_for_bam_samples.txt",
          ch_mismatching_metadata_sample_names.collect()
        )


        ch_bam_sample_qc_pass_matches_search.matched.map {
          it -> log.warn """Sample ${it} was provided with bam file. It was processed by pipeline before and was marked as QC_PASS by artic-ncov2019.
              Sample will be processed by pipeline again and will overwrite previous results
          """
        }
        store_notification_processed_already(
          "${workflow.start}-already_processed_bam_samples.txt",
          ch_bam_sample_qc_pass_matches_search.matched.collect()
        )

        ch_bam_matching_metadata.ifEmpty {
            log.error """\
              ERROR: No illumina bam file found matching samples, provided by I-SEHA sample metdata import.
                This may be caused by failure in loading sample metadata from .tsv file to the database, or metadata .tsv file not matching any bam files provided.
                Aborting!
            """
            System.exit(1)
        }

    emit:
        ch_bam_matching_metadata
}

// Process which acts as a workaround to append boolean value, whether the sample exists in another channel
process append_metadata_match_to_sample_file{
  input:
    tuple val(sample_name), file(bam)
    val samples_with_meta

  output:
    tuple val(sample_name), file(bam), val(has_meta)

  script:
    has_meta = samples_with_meta.contains(sample_name)

  """
  """
}

process store_bams_not_matched{
    publishDir COVID_PIPELINE_MISSING_METADATA_PATH, mode: 'copy', overwrite: true

    input:
      file(bam)

    output:
      file(bam)

    script:

    """
    """
}

