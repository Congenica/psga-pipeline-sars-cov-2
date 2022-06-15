/* This workflow runs per sample analysis run */

/*
 * Submit ncov sample results
 */
process submit_ncov_results {
  publishDir "${PSGA_OUTPUT_PATH}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: '*.{fa,png,csv}'

  input:
    path ncov_all_sample_results

  output:
    path "${output_path}", emit: ch_ncov_qc_all_samples_csv

  script:
    output_path = "ncov_qc.csv"

  """
  # Merge ncov QC results into one single file

  # extract the header from a lineage report
  sample_qc="\$(ls *.qc.csv | head -n 1)"
  # copy over the header only (first line)
  awk 'FNR == 1' "\${sample_qc}" > ${output_path}
  # copy over the record only from all files (second line)
  awk 'FNR == 2' *.qc.csv >> ${output_path}

  # remove the qc.csv files containing single samples
  rm *.qc.csv
  """
}


/*
 * Submit pangolin sample results
 */
process submit_pangolin_results {
  publishDir "${PSGA_OUTPUT_PATH}/pangolin", mode: 'copy', overwrite: true, pattern: 'all_lineages_report.csv'

  input:
    file input_dir

  output:
    path "${output_path}", emit: ch_pangolin_all_lineages

  script:
    output_path = "all_lineages_report.csv"

  """
  # Merge pangolin lineages into one single file

  # extract the header from a lineage report
  a_lineage_report="\$(ls *_lineage_report.csv | head -n 1)"
  # copy over the header only (first line)
  awk 'FNR == 1' "\${a_lineage_report}" > ${output_path}
  # copy over the record only from all files (second line)
  awk 'FNR == 2' *_lineage_report.csv >> ${output_path}
  """
}


/*
 * Concatenate ncov and pangolin CSV files
 */
process submit_pipeline_output_csv {
  publishDir "${PSGA_OUTPUT_PATH}/merged_output", mode: 'copy', overwrite: true, pattern: 'pipeline_output.csv'
  publishDir "${PSGA_OUTPUT_PATH}/notifications", mode: 'copy', overwrite: true, pattern: 'samples_{unknown,failed,passed}_{ncov_qc,pangolin}.txt'
  publishDir "${PSGA_OUTPUT_PATH}/logs", mode: 'copy', overwrite: true, pattern: '*.log'

  input:
    val ch_analysis_run_name
    path ch_metadata
    path ch_qc_ncov_result_csv_file
    path ch_pangolin_csv_file

  output:
    path ch_merged_csv_file, emit: ch_merged_csv_file
    path "*.txt", emit: ch_samples_files_by_qc
    path "*.log"

  script:
    ch_merged_csv_file = "pipeline_output.csv"

  """
  ncov_opt=""
  if [[ -f "${ch_qc_ncov_result_csv_file}" ]]; then
      ncov_opt="--ncov-qc-csv-file \"${ch_qc_ncov_result_csv_file}\""
  fi

  python ${PSGA_ROOT_PATH}/scripts/sars_cov_2/merge_ncov_pangolin_csv_files.py \
    --analysis-run-name "${ch_analysis_run_name}" \
    --metadata-file "${ch_metadata}" \
    --pangolin-csv-file "${ch_pangolin_csv_file}" \
    --merged-output-csv-file "${ch_merged_csv_file}" \
    \${ncov_opt}
  """
}


/*
 * Prepare and save the results for this analysis run.
 * This includes the output from ncov and pangolin
 */
workflow submit_analysis_run_results {
    take:
        ch_metadata
        ch_ncov_all_samples_results
        ch_pangolin_csvs
    main:

        if ( params.filetype == "fasta" ) {
            // mock ncov
            ch_ncov_submitted = ch_ncov_all_samples_results

        } else {
            // this will be bam or fastq

            submit_ncov_results(ch_ncov_all_samples_results.collect())

            ch_ncov_submitted = submit_ncov_results.out.ch_ncov_qc_all_samples_csv

        }

        submit_pangolin_results(ch_pangolin_csvs.collect())

        submit_pipeline_output_csv(
            params.run,
            ch_metadata,
            ch_ncov_submitted,
            submit_pangolin_results.out.ch_pangolin_all_lineages,
        )

        ch_analysis_run_results_submitted = submit_pipeline_output_csv.out.ch_merged_csv_file

    emit:
        ch_analysis_run_results_submitted
}
