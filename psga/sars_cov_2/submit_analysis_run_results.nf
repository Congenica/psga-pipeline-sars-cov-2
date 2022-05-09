/* This workflow runs per sample analysis run */

/*
 * Merge ncov QC results into one single file
 */
process merge_ncov2019_artic_qc_sample_files {
  input:
    file input_dir

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
    path ch_ncov_fasta
    path qc_plot_image
    path qc_csv

  output:
    path ch_ncov_fasta
    path qc_plot_image
    path qc_csv

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
    file ch_qc_ncov_result_csv_file
    file ch_qc_plot_files
    file ch_pangolin_all_lineages

  output:
    path ch_load_results_to_db_done

  script:
    directory_with_qc_depth_files = "./"
    ch_load_results_to_db_done = "load_results_to_db.done"

  """
  python ${PSGA_ROOT_PATH}/scripts/load_ncov_data_to_db.py \
    --ncov-qc-csv-file "${ch_qc_ncov_result_csv_file}" \
    --ncov-qc-depth-directory "${directory_with_qc_depth_files}" \
    --analysis-run-name "${ch_analysis_run_name}"

  python ${PSGA_ROOT_PATH}/scripts/load_pangolin_data_to_db.py \
    --pangolin-lineage-report-file "${ch_pangolin_all_lineages}" \
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
        analysis_run
        ch_ncov_qc_csvs
        ch_ncov_fastas
        ch_ncov_qc_plots
        ch_pangolin_csvs
    main:

        merge_ncov2019_artic_qc_sample_files(ch_ncov_qc_csvs)

        merge_pangolin_sample_files(ch_pangolin_csvs)

        store_ncov2019_artic_output(
            ch_ncov_fastas,
            ch_ncov_qc_plots,
            merge_ncov2019_artic_qc_sample_files.out.ch_ncov_qc_all_samples
        )

        store_pangolin_output(
            merge_pangolin_sample_files.out.ch_pangolin_all_lineages
        )

        ch_analysis_run_results_submitted = load_results_to_db(
            analysis_run,
            merge_ncov2019_artic_qc_sample_files.out.ch_ncov_qc_all_samples,
            ch_ncov_qc_plots,
            merge_pangolin_sample_files.out.ch_pangolin_all_lineages
        )

    emit:
        ch_analysis_run_results_submitted
}
