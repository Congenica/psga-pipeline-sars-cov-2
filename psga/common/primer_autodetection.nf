/*
 * Run: primer-autodetection (use first read only for illumina)
 */
process primer_autodetection {
  publishDir "${params.output_path}/primer_autodetection", mode: 'copy', overwrite: true, pattern: '{*_primer_data.csv,*_primer_detection.csv,*_trimming_report.txt}'

  tag "${task.index} - ${fastq}"

  input:
    path fastq
    val pathogen

  output:
    tuple path("*_primer.txt"), path(fastq), emit: ch_files
    path "*_primer_data.csv", emit: ch_primer_data
    path "*_primer_detection.csv", emit: ch_primer_coverage
    path "*_trimming_report.txt", emit: ch_trimming_report

  shell:
  '''
  # select the input fastq file based on illumina or ont and trim adaptor sequences.
  # For illumina, 1 read is sufficient for primer autodetection
  if [[ "!{params.sequencing_technology}" == "illumina" ]]; then
    file_1="$(ls *_1.fastq.gz)"
    file_2="$(ls *_2.fastq.gz)"
    sample_id=$( echo ${file_1} | cut -d '_' -f1)
    trimmed_file_1="${sample_id}_1_val_1.fq.gz"
    trim_galore --paired --gzip ${file_1} ${file_2} > ${sample_id}_trimming_report.txt 2>&1
  elif [[ "!{params.sequencing_technology}" == "ont" ]]; then
    file_1="!{fastq}"
    sample_id=$( echo ${file_1} | cut -d '_' -f1)
    trimmed_file_1="${sample_id}_trimmed.fq.gz"
    porechop --extra_end_trim !{params.porechop_extra_end_trim} --format fastq.gz -i ${file_1} -o ${trimmed_file_1} > ${sample_id}_trimming_report.txt 2>&1
  else
    echo "sequencing technology must be either 'illumina' or 'ont' for primer autodetection"
    exit 1
  fi

  # extract the sample id
  primer_index="/primer_schemes/!{pathogen}_primer_index.csv"

  python ${PSGA_ROOT_PATH}/scripts/common/primer_autodetection.py \
    --primer-index "${primer_index}" \
    --sample-fastq "${trimmed_file_1}" \
    --sample-id "${sample_id}" \
    --primer-input !{params.kit} \
    --primer-base-window !{params.primer_base_window}
  '''
}
