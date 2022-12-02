/* This workflow runs per sample analysis run */

include { submit_results as submit_contamination_removal_results } from './common/submit_results.nf'
include { submit_results as submit_primer_autodetection_results } from './common/submit_results.nf'
include { submit_results as submit_ncov_qc_results } from './common/submit_results.nf'
include { submit_results as submit_pangolin_results } from './common/submit_results.nf'


/*
 * Submit the pipeline output files
 */
process submit_pipeline_results_files {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true, pattern: 'result{s.csv,s.json,files.json}'
  publishDir "${params.output_path}/logs", mode: 'copy', overwrite: true, pattern: '*.log'

  input:
    path ch_metadata
    path ch_contamination_removal_csv_file
    path ch_primer_autodetection_csv_file
    path ch_ncov_qc_result_csv_file
    path ch_pangolin_csv_file

  output:
    path "results.csv", emit: ch_results_csv
    path "results.json", emit: ch_results_json
    path "resultfiles.json", emit: ch_resultfiles_json
    path "*.log"

  script:
  """
  ncov_opt=""
  if [[ -f "${ch_contamination_removal_csv_file}" ]] \
     && [[ -f "${ch_primer_autodetection_csv_file}" ]] \
     && [[ -f "${ch_ncov_qc_result_csv_file}" ]]; then
      ncov_opt="--contamination-removal-csv-file \"${ch_contamination_removal_csv_file}\" --primer-autodetection-csv-file \"${ch_primer_autodetection_csv_file}\" --ncov-qc-csv-file \"${ch_ncov_qc_result_csv_file}\""
  fi

  python ${PSGA_ROOT_PATH}/scripts/sars_cov_2/generate_pipeline_results_files.py \
    --analysis-run-name "${params.run}" \
    --metadata-file "${ch_metadata}" \
    --pangolin-csv-file "${ch_pangolin_csv_file}" \
    --output-path "${params.output_path}" \
    --sequencing-technology "${params.sequencing_technology}" \
    \${ncov_opt}
  """
}


/*
 * Prepare and save the results for this analysis run.
 */
workflow submit_analysis_run_results {
    take:
        ch_metadata
        ch_contamination_removal_csvs
        ch_primer_autodetection_csvs
        ch_ncov_qc_csvs
        ch_pangolin_csvs
    main:

        if ( params.sequencing_technology == "unknown" ) {
            // ncov was not executed. Therefore, mock contamination removal, primer_autodetection and ncov
            ch_contamination_removal_submitted = ch_contamination_removal_csvs
            ch_primer_autodetection_submitted = ch_primer_autodetection_csvs
            ch_ncov_qc_submitted = ch_ncov_qc_csvs
        } else {
            ch_contamination_removal_submitted = submit_contamination_removal_results(
                ch_contamination_removal_csvs,
                'contamination_removal.csv',
                'sample_id',
                'contamination_removal'
            )
            ch_primer_autodetection_submitted = submit_primer_autodetection_results(
                ch_primer_autodetection_csvs,
                'primer_autodetection.csv',
                'sample_id',
                'primer_autodetection'
            )
            ch_ncov_qc_submitted = submit_ncov_qc_results(
                ch_ncov_qc_csvs,
                'ncov_qc.csv',
                'sample_name',
                'ncov2019-artic'
            )
        }

        ch_pangolin_submitted = submit_pangolin_results(
            ch_pangolin_csvs,
            'all_lineages_report.csv',
            'taxon',
            'pangolin'
        )

        // if the channel is empty, ifEmpty(file(<default_header>)) is used
        // so that this process is always triggered
        submit_pipeline_results_files(
            ch_metadata,
            ch_contamination_removal_submitted.ifEmpty(file(params.contamination_removal_empty_csv)),
            ch_primer_autodetection_submitted.ifEmpty(file(params.primer_autodetection_empty_csv)),
            ch_ncov_qc_submitted.ifEmpty(file(params.ncov_qc_empty_csv)),
            ch_pangolin_submitted.ifEmpty(file(params.pangolin_empty_csv)),
        )

        ch_analysis_run_results_submitted = submit_pipeline_results_files.out.ch_results_json

    emit:
        ch_analysis_run_results_submitted
}
