/*
 * Generate synthetic output
 */
process generate_synthetic_output {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true, pattern: 'result{s.csv,files.json}'

  input:
    path ch_metadata
    path input_files

  output:
    path ch_output_csv_file, emit: ch_output_csv_file
    path ch_output_json_file, emit: ch_output_json_file

  script:
    ch_output_csv_file = "results.csv"
    ch_output_json_file = "resultfiles.json"

  """
  output_csv_file="results.csv"
  output_json_file="resultfiles.json"

  python ${PSGA_ROOT_PATH}/scripts/synthetic/generate_pipeline_results_files.py \
    --metadata-file "${ch_metadata}" \
    --sequencing-technology "${params.sequencing_technology}" \
    --output-path "${params.output_path}"
  """
}
