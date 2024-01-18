/*
 * Generic process for concatenating and submitting sample results.
 */
process submit_results {
  publishDir "${params.output_path}/${ch_publish_dirname}", mode: 'copy', overwrite: true

  input:
    path ch_input_csvs
    val ch_output_path
    val ch_sortby_col
    val ch_publish_dirname

  output:
    path "${ch_output_path}", emit: ch_all_samples_csv

  script:
  """
  python ${PSGA_ROOT_PATH}/scripts/concat_csv.py --input-path . --output-csv-path ${ch_output_path} --sortby-col ${ch_sortby_col}
  """
}
