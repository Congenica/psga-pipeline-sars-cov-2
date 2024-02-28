/*
 * Run: read-it-and-keep
 * https://github.com/GlobalPathogenAnalysisService/read-it-and-keep
 * This tool keeps the reads that match the provided target genome.
 */

process CONTAMINATION_REMOVAL {
  publishDir "${params.output_path}/contamination_removal", mode: 'copy', overwrite: true, pattern: '{*_contamination_removal.csv,cleaned_fastq/*.fastq.gz,counting/*.txt}'

  input:
    val ref_genome_fasta
    tuple val(meta), path(read_paths)

  output:
    path "*_contamination_removal.csv", emit: ch_contamination_removal_csv
    path "cleaned_fastq/*.fastq.gz", emit: ch_output_file
    path "counting/*.txt"

    tuple val(meta), path("cleaned_fastq/${sample_id}_*.fastq.gz"), emit: ch_cleaned_fastq

  script:

    // Common Variables to both
    sample_id = meta.SAMPLE_ID
    rik_output_file = "counting/${sample_id}.txt"
    // Note: Added _1 to ont file here, it shouldn't matter...
    // As we will use the path and not reconstruct it elsewhere
    cleaned_fastq_file_1 = "cleaned_fastq/${sample_id}_1.fastq.gz"
    output_csv = "${sample_id}_contamination_removal.csv"

    // NOTE: readItAndKeep always compresses the output, not matter what the input was so assume filename is .gz
    // This was gzipping the input file, but it doesn't seem to be necessary
    if( params.sequencing_technology == 'illumina' ) {
      // Only if illumina
      // TODO: ENSURE Ordering is correctly ordered
      reads_file_1 = read_paths[0]
      reads_file_2 = read_paths[1]
      cleaned_fastq_file_2 = "cleaned_fastq/${sample_id}_2.fastq.gz"
      """
      mkdir -p cleaned_fastq counting

      readItAndKeep \
        --tech ${params.sequencing_technology} \
        --ref_fasta ${ref_genome_fasta} \
        --reads1 ${reads_file_1} \
        --reads2 ${reads_file_2} \
        --outprefix out 2>&1 \
      | tee ${rik_output_file}

      mv -f out.reads_1.fastq.gz ${cleaned_fastq_file_1}
      mv -f out.reads_2.fastq.gz ${cleaned_fastq_file_2}

      python ${PSGA_ROOT_PATH}/scripts/contamination_removal.py \
        --input-path "${rik_output_file}" \
        --output-csv-path "${output_csv}" \
        --sample-id "${sample_id}"
      """
    }

    else if( params.sequencing_technology == 'ont' ) {
      reads_file_1 = read_paths
      """
      mkdir -p cleaned_fastq counting

      readItAndKeep \
        --tech ${params.sequencing_technology} \
        --ref_fasta ${ref_genome_fasta} \
        --reads1 ${reads_file_1} \
        --outprefix out 2>&1 \
      | tee ${rik_output_file}

      mv -f out.reads.fastq.gz ${cleaned_fastq_file_1}

      python ${PSGA_ROOT_PATH}/scripts/contamination_removal.py \
        --input-path "${rik_output_file}" \
        --output-csv-path "${output_csv}" \
        --sample-id "${sample_id}"
      """
    }
}