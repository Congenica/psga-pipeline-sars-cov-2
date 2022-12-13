/*
 * Generate s_aureus output
 */
process generate_s_aureus_output {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true, pattern: 'result{s.csv,files.json}'

  input:
    path ch_metadata
    path input_files
    path annotation_summary
    path variants_txt_for_csv_file
    path all_software_versions
    path checkm_results_txt
    path assembly_json, stageAs: 'assembly.json'
    path antimicrobial_protein_report
    path antimicrobial_gene_report
    path mykrobe_versions, stageAs: 'mykrobe_versions.yml'
    path mlst_versions, stageAs: 'mlst_versions.yml'

  output:
    path ch_output_csv_file, emit: ch_output_csv_file
    path ch_output_json_file, emit: ch_output_json_file

  script:
    ch_output_csv_file = "results.csv"
    ch_output_json_file = "resultfiles.json"

  """
  output_csv_file="results.csv"
  output_json_file="resultfiles.json"

  python ${PSGA_ROOT_PATH}/scripts/s_aureus/generate_results.py \
    --metadata-file "${ch_metadata}" \
    --output-csv-file "${ch_output_csv_file}" \
    --output-json-file "${ch_output_json_file}" \
    --output-path "${params.output_path}" \
    --sequencing-technology "${params.sequencing_technology}"
  """
}
