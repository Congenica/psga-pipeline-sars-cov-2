/*
 * Run: read-it-and-keep
 * https://github.com/GlobalPathogenAnalysisService/read-it-and-keep
 * This tool keeps the reads that match the provided target genome.
 */
process contamination_removal_ont {
  publishDir "${params.output_path}/contamination_removal", mode: 'copy', overwrite: true, pattern: '{*_contamination_removal.csv,cleaned_fastq/*.fastq.gz,counting/*.txt}'

  tag "${task.index} - ${file_1}"

  input:
    val ref_genome_fasta
    path file_1

  output:
    path "*_contamination_removal.csv", emit: ch_contamination_removal_csv
    path "cleaned_fastq/*.fastq.gz", emit: ch_output_file
    path "counting/*.txt"

  shell:
  '''
  # NOTE: readItAndKeep always compresses the output, not matter what the input was
  # therefore, make sure that the input file is compressed, in order to have the right output file name

  ref_genome_fasta=!{ref_genome_fasta}
  file_1=!{file_1}
  sample_id=$( echo ${file_1} | cut -d '_' -f1 )
  rik_output_file="counting/${sample_id}.txt"
  output_csv="${sample_id}_contamination_removal.csv"

  if [ "${file_1##*.}" != "gz" ]; then
      # keep input file, otherwise nextflow complains
      gzip -c ${file_1} > ${file_1}.gz
      file_1="${file_1}.gz"
  fi

  mkdir -p cleaned_fastq counting

  readItAndKeep --tech ont --ref_fasta ${ref_genome_fasta} --reads1 ${file_1} --outprefix out 2>&1 | tee ${rik_output_file}

  mv -f out.reads.fastq.gz cleaned_fastq/${sample_id}_1.fastq.gz

  python ${PSGA_ROOT_PATH}/scripts/contamination_removal.py \
    --input-path "${rik_output_file}" \
    --output-csv-path "${output_csv}" \
    --sample-id "${sample_id}"
  '''
}

process contamination_removal_illumina {
  publishDir "${params.output_path}/contamination_removal", mode: 'copy', overwrite: true, pattern: '{*_contamination_removal.csv,cleaned_fastq/*.fastq.gz,counting/*.txt}'

  tag "${task.index} - [${file_1}, ${file_2}]"

  input:
    val ref_genome_fasta
    tuple path(file_1), path(file_2)

  output:
    path "*_contamination_removal.csv", emit: ch_contamination_removal_csv
    path "cleaned_fastq/*.fastq.gz", emit: ch_output_file
    path "counting/*.txt"

  shell:
  '''
  # NOTE: readItAndKeep always compresses the output, not matter what the input was
  # therefore, make sure that the input file is compressed, in order to have the right output file name

  ref_genome_fasta=!{ref_genome_fasta}
  file_1=!{file_1}
  file_2=!{file_2}
  sample_id=$( echo ${file_1} | cut -d '_' -f1 )
  rik_output_file="counting/${sample_id}.txt"
  output_csv="${sample_id}_contamination_removal.csv"

  if [ "${file_1##*.}" != "gz" ]; then
      # keep input file, otherwise nextflow complains
      gzip -c ${file_1} > ${file_1}.gz
      file_1="${file_1}.gz"
  fi
  if [ "${file_2##*.}" != "gz" ]; then
      # keep input file, otherwise nextflow complains
      gzip -c ${file_2} > ${file_2}.gz
      file_2="${file_2}.gz"
  fi

  mkdir -p cleaned_fastq counting

  readItAndKeep --tech illumina --ref_fasta ${ref_genome_fasta} --reads1 ${file_1} --reads2 ${file_2} --outprefix out 2>&1 | tee ${rik_output_file}

  mv -f out.reads_1.fastq.gz cleaned_fastq/${sample_id}_1.fastq.gz
  mv -f out.reads_2.fastq.gz cleaned_fastq/${sample_id}_2.fastq.gz

  python ${PSGA_ROOT_PATH}/scripts/contamination_removal.py \
    --input-path "${rik_output_file}" \
    --output-csv-path "${output_csv}" \
    --sample-id "${sample_id}"
  '''
}
