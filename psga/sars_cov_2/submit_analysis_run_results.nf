/* This workflow runs per sample analysis run */

/*
 * Submit ncov sample results
 */
process submit_ncov_results {
  publishDir "${params.output_path}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: 'ncov_qc.csv'

  input:
    path ch_ncov_qc_csvs

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
  """
}


/*
 * Submit pangolin sample results
 */
process submit_pangolin_results {
  publishDir "${params.output_path}/pangolin", mode: 'copy', overwrite: true, pattern: 'all_lineages_report.csv'

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
 * Submit the merged ncov and pangolin output file
 */
process submit_pipeline_results_files {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true, pattern: 'result{s.csv,files.json}'
  publishDir "${params.output_path}/notifications", mode: 'copy', overwrite: true, pattern: 'samples_{unknown,failed,passed}_{ncov_qc,pangolin}.txt'
  publishDir "${params.output_path}/logs", mode: 'copy', overwrite: true, pattern: '*.log'

  input:
    val ch_analysis_run_name
    path ch_metadata
    path ch_qc_ncov_result_csv_file
    path ch_pangolin_csv_file

  output:
    path ch_output_csv_file, emit: ch_output_csv_file
    path ch_output_json_file, emit: ch_output_json_file
    path "*.txt", emit: ch_samples_files_by_qc
    path "*.log"

  script:
    ch_output_csv_file = "results.csv"
    ch_output_json_file = "resultfiles.json"

  """
  ncov_opt=""
  if [[ -f "${ch_qc_ncov_result_csv_file}" ]]; then
      ncov_opt="--ncov-qc-csv-file \"${ch_qc_ncov_result_csv_file}\""
  fi

  python ${PSGA_ROOT_PATH}/scripts/sars_cov_2/generate_pipeline_results_files.py \
    --analysis-run-name "${ch_analysis_run_name}" \
    --metadata-file "${ch_metadata}" \
    --pangolin-csv-file "${ch_pangolin_csv_file}" \
    --output-csv-file "${ch_output_csv_file}" \
    --output-json-file "${ch_output_json_file}" \
    --output-path "${params.output_path}" \
    --sequencing-technology "${params.sequencing_technology}" \
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
        ch_ncov_qc_csvs
        ch_pangolin_csvs
    main:

        if ( params.sequencing_technology == "unknown" ) {
            // ncov was not executed. Therefore, mock ncov
            ch_ncov_submitted = ch_ncov_qc_csvs
        } else {
            submit_ncov_results(ch_ncov_qc_csvs)
            ch_ncov_submitted = submit_ncov_results.out.ch_ncov_qc_all_samples_csv
        }

        submit_pangolin_results(ch_pangolin_csvs.collect())

        // if the channel is empty, ifEmpty(file(<default_header>)) is used
        // so that this process is always triggered
        submit_pipeline_results_files(
            params.run,
            ch_metadata,
            ch_ncov_submitted.ifEmpty(file(params.ncov_qc_empty_csv)),
            submit_pangolin_results.out.ch_pangolin_all_lineages.ifEmpty(file(params.pangolin_empty_csv)),
        )

        ch_analysis_run_results_submitted = submit_pipeline_results_files.out.ch_output_csv_file

    emit:
        ch_analysis_run_results_submitted
}
