process PRIMER_AUTODETECTION {
  publishDir "${params.output_path}/primer_autodetection", mode: 'copy', overwrite: true, pattern: '{*_primer_data.csv,*_primer_detection.csv}'

  input:
    tuple val(meta), path(read_paths)

  output:
    path "*_primer_data.csv", emit: ch_primer_data
    path "*_primer_detection.csv", emit: ch_primer_coverage

    // Primer is written to a file so it can be added to metadata
    tuple val(meta), path(read_paths), path("${sample_id}_primer.txt"), emit: ch_primer_detected

  script:

    sample_id = meta.SAMPLE_ID
    primer_index="/primer_schemes/SARS-CoV-2_primer_index.csv"

    if (params.sequencing_technology == "illumina") {
      read_path = read_paths[0]
    } else {
      read_path = read_paths
    }

    """
    python ${PSGA_ROOT_PATH}/scripts/primer_autodetection.py \
      --primer-index "${primer_index}" \
      --sample-fastq "${read_path}" \
      --sample-id "${sample_id}" \
      --primer-input ${params.kit}
    """
}
