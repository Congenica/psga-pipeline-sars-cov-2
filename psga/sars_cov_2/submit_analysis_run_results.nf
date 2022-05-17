/* This workflow runs per sample analysis run */

include { store_notification as store_missing_ncov_qc_notification } from '../common/utils.nf'
include { store_notification as store_failed_ncov_qc_notification } from '../common/utils.nf'
include { store_notification as store_passed_ncov_qc_notification } from '../common/utils.nf'
include { store_notification as store_unknown_pangolin_notification } from '../common/utils.nf'
include { store_notification as store_failed_pangolin_notification } from '../common/utils.nf'
include { store_notification as store_passed_pangolin_notification } from '../common/utils.nf'

/*
 * Merge ncov QC results into one single file
 */
process merge_ncov2019_artic_qc_sample_files {
  input:
    path ncov_all_sample_results

  output:
    path "${output_path}", emit: ch_ncov_qc_all_samples

  script:
    output_path = "ncov_qc.csv"

  """
  # extract the header from a lineage report
  sample_qc="\$(ls *.qc.csv | head -n 1)"
  # copy over the header only (first line)
  awk 'FNR == 1' "\${sample_qc}" > ${output_path}
  # copy over the record only from all files (second line)
  awk 'FNR == 2' *.qc.csv >> ${output_path}
  """
}

/*
 * Store ncov2019_artic output
 */
process store_ncov2019_artic_output {
  publishDir "${PSGA_OUTPUT_PATH}/ncov2019-artic", mode: 'copy', overwrite: true

  input:
    path ncov_all_sample_results
    path all_sample_qc_csv

  output:
    path ncov_all_sample_results
    path all_sample_qc_csv

  script:

  """
  """
}

/*
 * Merge pangolin lineages into one single file
 */
process merge_pangolin_sample_files {
  input:
    file input_dir

  output:
    path "${output_path}", emit: ch_pangolin_all_lineages

  script:
    output_path = "all_lineages_report.csv"

  """
  # extract the header from a lineage report
  a_lineage_report="\$(ls *_lineage_report.csv | head -n 1)"
  # copy over the header only (first line)
  awk 'FNR == 1' "\${a_lineage_report}" > ${output_path}
  # copy over the record only from all files (second line)
  awk 'FNR == 2' *_lineage_report.csv >> ${output_path}
  """
}

/*
 * Store pangolin output
 */
process store_pangolin_output {
  publishDir "${PSGA_OUTPUT_PATH}/pangolin", mode: 'copy', overwrite: true

  input:
    path ch_pangolin_csv

  output:
    path ch_pangolin_csv

  script:

  """
  """
}

/*
 * Load sample results to the database.
 * Use one single process to enforce a sequential loading.
 */
process load_results_to_db {
  input:
    val ch_analysis_run_name
    path ncov_all_sample_results
    path ch_qc_ncov_result_csv_file
    path ch_pangolin_all_lineages

  output:
    path ch_load_results_to_db_done, emit: ch_load_results_to_db_done
    path ch_samples_without_ncov_qc, emit: ch_samples_without_ncov_qc
    path ch_samples_with_failed_ncov_qc, emit: ch_samples_with_failed_ncov_qc
    path ch_samples_with_passed_ncov_qc, emit: ch_samples_with_passed_ncov_qc
    path ch_samples_with_unknown_pangolin_status, emit: ch_samples_with_unknown_pangolin_status
    path ch_samples_with_failed_pangolin_status, emit: ch_samples_with_failed_pangolin_status
    path ch_samples_with_passed_pangolin_status, emit: ch_samples_with_passed_pangolin_status

  script:
    directory_with_qc_depth_files = "./"
    ch_load_results_to_db_done = "load_results_to_db.done"
    ch_samples_without_ncov_qc = "samples_without_ncov_qc.txt"
    ch_samples_with_failed_ncov_qc = "samples_with_failed_ncov_qc.txt"
    ch_samples_with_passed_ncov_qc = "samples_with_passed_ncov_qc.txt"
    ch_samples_with_unknown_pangolin_status = "samples_with_unknown_pangolin_status.txt"
    ch_samples_with_failed_pangolin_status = "samples_with_failed_pangolin_status.txt"
    ch_samples_with_passed_pangolin_status = "samples_with_passed_pangolin_status.txt"

  """
  # Note: the samples in load_pangolin could be a subset of the samples in load_ncov.
  # This can happen if ncov QC failed or ncov failed to run for a sample in general.

  python ${PSGA_ROOT_PATH}/scripts/load_ncov_data_to_db.py \
    --ncov-qc-csv-file "${ch_qc_ncov_result_csv_file}" \
    --ncov-qc-depth-directory "${directory_with_qc_depth_files}" \
    --samples-without-ncov-qc-file "${ch_samples_without_ncov_qc}" \
    --samples-with-failed-ncov-qc-file "${ch_samples_with_failed_ncov_qc}" \
    --samples-with-passed-ncov-qc-file "${ch_samples_with_passed_ncov_qc}" \
    --analysis-run-name "${ch_analysis_run_name}"

  python ${PSGA_ROOT_PATH}/scripts/load_pangolin_data_to_db.py \
    --pangolin-lineage-report-file "${ch_pangolin_all_lineages}" \
    --samples-with-unknown-pangolin-status-file "${ch_samples_with_unknown_pangolin_status}" \
    --samples-with-failed-pangolin-status-file "${ch_samples_with_failed_pangolin_status}" \
    --samples-with-passed-pangolin-status-file "${ch_samples_with_passed_pangolin_status}" \
    --analysis-run-name "${ch_analysis_run_name}"

  touch ${ch_load_results_to_db_done}
  """
}


/*
 * Prepare and save the results for this analysis run.
 * This includes the output from ncov and pangolin
 */
workflow submit_analysis_run_results {
    take:
        ch_ncov_all_samples_results
        ch_pangolin_csvs
    main:

        // split by qc_csv so that only the relevant files are passed to the downstream processes
        ch_ncov_all_samples_results
            .flatten()
            .branch {
                qc_csv_only: it =~ /^.*\.qc\.csv$/
                    return it
                filter_out_qc_csv: true
                    return it
            }
            .set{ ch_ncov_all_samples_results_branch }

        merge_ncov2019_artic_qc_sample_files(ch_ncov_all_samples_results_branch.qc_csv_only.collect())

        merge_pangolin_sample_files(ch_pangolin_csvs)

        store_ncov2019_artic_output(
            ch_ncov_all_samples_results_branch.filter_out_qc_csv.collect(),
            merge_ncov2019_artic_qc_sample_files.out.ch_ncov_qc_all_samples
        )

        store_pangolin_output(
            merge_pangolin_sample_files.out.ch_pangolin_all_lineages
        )

        load_results_to_db(
            params.run,
            ch_ncov_all_samples_results_branch.filter_out_qc_csv.collect(),
            merge_ncov2019_artic_qc_sample_files.out.ch_ncov_qc_all_samples,
            merge_pangolin_sample_files.out.ch_pangolin_all_lineages
        )

        ch_analysis_run_results_submitted = load_results_to_db.out.ch_load_results_to_db_done

        ch_missing_ncov_notification_submitted = store_missing_ncov_qc_notification(
            load_results_to_db.out.ch_samples_without_ncov_qc
        )
        ch_failed_ncov_notification_submitted = store_failed_ncov_qc_notification(
            load_results_to_db.out.ch_samples_with_failed_ncov_qc
        )
        ch_passed_ncov_notification_submitted = store_passed_ncov_qc_notification(
            load_results_to_db.out.ch_samples_with_passed_ncov_qc
        )
        ch_unknown_pangolin_notification_submitted = store_unknown_pangolin_notification(
            load_results_to_db.out.ch_samples_with_unknown_pangolin_status
        )
        ch_failed_pangolin_notification_submitted = store_failed_pangolin_notification(
            load_results_to_db.out.ch_samples_with_failed_pangolin_status
        )
        ch_passed_pangolin_notification_submitted = store_passed_pangolin_notification(
            load_results_to_db.out.ch_samples_with_passed_pangolin_status
        )

    emit:
        ch_analysis_run_results_submitted
        ch_missing_ncov_notification_submitted
        ch_failed_ncov_notification_submitted
        ch_passed_ncov_notification_submitted
        ch_unknown_pangolin_notification_submitted
        ch_failed_pangolin_notification_submitted
        ch_passed_pangolin_notification_submitted
}