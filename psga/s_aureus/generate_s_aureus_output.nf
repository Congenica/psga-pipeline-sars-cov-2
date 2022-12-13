/*
 * Generate s_aureus output
 */
process generate_s_aureus_output {
  publishDir "${params.output_path}", mode: 'copy', overwrite: true

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
    path genome_assembly_name

  output:
    path "*-results.csv", emit: ch_output_csv_file
    path "*-resultfiles.json", emit: ch_output_json_file

  script:
  """
  python ${PSGA_ROOT_PATH}/scripts/s_aureus/generate_results.py \
    --output-path "${params.output_path}"
  """
}
