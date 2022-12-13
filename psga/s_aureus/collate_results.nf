/*
 * Combine all of the results.csv files and resultfiles.json
 */
process collate_results {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true

  input:
    path ch_output_csv_file
    path ch_output_json_file

  output:
    path "results.csv", emit: global_csv_file
    path "resultfiles.jaon", emit: global_json_file

  script:
  """
    python ${PSGA_ROOT_PATH}/scripts/s_aureus/collate_results.py
  """
}