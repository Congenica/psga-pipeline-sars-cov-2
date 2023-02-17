/*
 * Run: primer-autodetection (use first read only for illumina)
 */
process primer_autodetection {
  publishDir "${params.output_path}/primer_autodetection", mode: 'copy', overwrite: true, pattern: '{*_primer_data.csv,*_primer_detection.csv}'

  tag "${task.index} - ${fastq}"

  input:
    path fastq
    val pathogen

  output:
    tuple path("*_primer.txt"), path(fastq), emit: ch_files
    path "*_primer_data.csv", emit: ch_primer_data
    path "*_primer_detection.csv", emit: ch_primer_coverage

  shell:
  '''
  # select the input fastq file based on illumina or ont.
  # For illumina, 1 read is sufficient for primer autodetection
  if [[ "!{params.sequencing_technology}" == "illumina" ]]; then
    file_1="$(ls *_1.fastq.gz)"
  elif [[ "!{params.sequencing_technology}" == "ont" ]]; then
    file_1="!{fastq}"
  else
    echo "sequencing technology must be either 'illumina' or 'ont' for primer autodetection"
    exit 1
  fi

  # extract the sample id, whereas this is ID.fastq.gz or ID_1.fastq.gz
  sample_id=$( echo ${file_1} | cut -d '.' -f1 | cut -d '_' -f1)
  primer_index="/primer_schemes/!{pathogen}_primer_index.csv"

  python ${PSGA_ROOT_PATH}/scripts/common/primer_autodetection.py --primer-index "${primer_index}" --sample-fastq "${file_1}" --sample-id "${sample_id}" --primer-input !{params.kit}
  '''
}
