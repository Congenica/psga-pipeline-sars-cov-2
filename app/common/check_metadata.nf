/*
 * Check the metadata in the database
 */
process check_metadata {
  publishDir "${params.output_path}/logs", mode: 'copy', overwrite: true, pattern: '*.log'

  tag "${metadata}"

  input:
    path metadata

  output:
    path metadata, emit: ch_metadata
    path "*.log"

  shell:
  '''
  python ${PSGA_ROOT_PATH}/scripts/check_metadata.py \
    --metadata-path "!{metadata}" \
    --analysis-run-name "!{params.run}" \
    --sequencing-technology "!{params.sequencing_technology}"
  '''
}
