process PRIMER_AUTODETECTION {
  publishDir "${params.output_path}/primer_autodetection", mode: 'copy', overwrite: true, pattern: '{*_primer_data.csv,*_primer_detection.csv}'

  input:
    tuple val(meta), path(reads_file)

  output:
    path "*_primer_data.csv", emit: ch_primer_data
    path "*_primer_detection.csv", emit: ch_primer_coverage

    // Primer is written to a file so it can be added to metadata
    tuple val(meta), path(reads_file), path("${sample_id}_primer.txt"), emit: ch_primer_detected

  script:

    sample_id = meta.SAMPLE_ID
    primer_index="/primer_schemes/SARS-CoV-2_primer_index.csv"

    """
    python ${PSGA_ROOT_PATH}/scripts/primer_autodetection.py \
      --primer-index "${primer_index}" \
      --sample-fastq "${reads_file}" \
      --sample-id "${sample_id}" \
      --primer-input ${params.kit}
    """
}
