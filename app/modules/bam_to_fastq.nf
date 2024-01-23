process BAM_TO_FASTQ_ONT {
  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path(fastq_path)

  script:
    fastq_path = "${meta.SAMPLE_ID}.fastq"
  """
  # write singleton reads to (decompressed) fastq file. Do not append /1 and 2/ to the read name
  samtools bam2fq -nO ${bam} > ${fastq_path}
  """
}
