/*
 * Run: read-it-and-keep
 * https://github.com/GlobalPathogenAnalysisService/read-it-and-keep
 * This tool keeps the reads that match the provided target genome.
 */
process PROCESS_ONT_FASTQ {
  // Contamination Removal
  publishDir "${params.output_path}/contamination_removal", mode: 'copy', overwrite: true, pattern: '{*_contamination_removal.csv,cleaned_fastq/*.fastq.gz,counting/*.txt}'
  // FASTQC
  publishDir "${params.output_path}/fastqc", mode: 'copy', overwrite: true, pattern: '*_fastqc.zip'

  input:
    val ref_genome_fasta
    tuple val(meta), val(read_paths)

  output:
    // Contamination Removal
    path "*_contamination_removal.csv", emit: ch_contamination_removal_csv
    path "cleaned_fastq/*.fastq.gz", emit: ch_output_file
    path "counting/*.txt"
    // Fastqc
    path "*_fastqc.zip", emit: ch_fastqc_zip_report
    // primer autodetection
    tuple path("*_primer.txt"), path(fastq), emit: ch_files
    path "*_primer_data.csv", emit: ch_primer_data
    path "*_primer_detection.csv", emit: ch_primer_coverage

    tuple val(meta), path(cleaned_fastq_file), emit: ch_cleaned_fastq

  script:

    sample_id = meta.SAMPLE_ID
    reads_file = read_paths.seq_file_1
    rik_output_file = "counting/${sample_id}.txt"
    output_csv = "${sample_id}_contamination_removal.csv"
    cleaned_fastq_file = "cleaned_fastq/${sample_id}_1.fastq.gz"

    // NOTE: readItAndKeep always compresses the output, not matter what the input was so assume filename is .gz
    // This was gzipping the input file, but it doesn't seem to be necessary

    """
    mkdir -p cleaned_fastq counting

    readItAndKeep --tech ont --ref_fasta ${ref_genome_fasta} --reads1 ${reads_file} --outprefix out 2>&1 | tee ${rik_output_file}

    mv -f out.reads.fastq.gz ${cleaned_fastq_file}

    python ${PSGA_ROOT_PATH}/scripts/contamination_removal.py \
      --input-path "${rik_output_file}" \
      --output-csv-path "${output_csv}" \
      --sample-id "${sample_id}"

    fastqc -q --limits /limits.txt --outdir \$PWD ${cleaned_fastq_file}

    primer_index="/primer_schemes/sars_cov_2_primer_index.csv"

    python ${PSGA_ROOT_PATH}/scripts/primer_autodetection.py \
      --primer-index "${primer_index}" \
      --sample-fastq "${cleaned_fastq_file}" \
      --sample-id "${sample_id}" \
      --primer-input ${params.kit}
    """
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
