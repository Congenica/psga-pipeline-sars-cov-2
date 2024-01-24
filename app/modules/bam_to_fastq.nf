process BAM_TO_FASTQ_ONT {
  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), val(reads_map)

  script:
    fastq_path = "${meta.SAMPLE_ID}.fastq"
    reads_map = ["seq_file_1": file(fastq_path)]
  """
  # write singleton reads to (decompressed) fastq file. Do not append /1 and 2/ to the read name
  samtools bam2fq -nO ${bam} > ${fastq_path}
  """
}
