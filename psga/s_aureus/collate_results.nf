/*
 * Combine all of the results.csv files and resultfiles.json
 */
process collate_results {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true

  input:
    path ch_output_csv_file
    path ch_output_json_file
    path ch_metadata

  output:
    path "results.csv", emit: global_csv_file
    path "resultfiles.json", emit: global_json_file

  script:
  """
    # The framework copies and renames the metadata file for all samples. These files are identical and we only need one
    python ${PSGA_ROOT_PATH}/scripts/s_aureus/collate_results.py --metadata-file ${ch_metadata}
  """
}
