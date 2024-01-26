/*
 * Run: read-it-and-keep
 * https://github.com/GlobalPathogenAnalysisService/read-it-and-keep
 * This tool keeps the reads that match the provided target genome.
 */
process CONTAMINATION_REMOVAL {
  publishDir "${params.output_path}/contamination_removal", mode: 'copy', overwrite: true, pattern: '{*_contamination_removal.csv,cleaned_fastq/*.fastq.gz,counting/*.txt}'

  input:
    val ref_genome_fasta
    tuple val(meta), val(read_paths)

  output:
    path "*_contamination_removal.csv", emit: ch_contamination_removal_csv
    path "cleaned_fastq/*.fastq.gz", emit: ch_output_file
    path "counting/*.txt"

    tuple val(meta), path(cleaned_fastq_file), emit: ch_cleaned_fastq

  script:

    sample_id = meta.SAMPLE_ID
    reads_file = read_paths.seq_file_1
    rik_output_file = "counting/${sample_id}.txt"
    cleaned_fastq_file = "cleaned_fastq/${sample_id}_1.fastq.gz"
    output_csv = "${sample_id}_contamination_removal.csv"

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
    """
}