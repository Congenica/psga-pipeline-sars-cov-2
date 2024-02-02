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


process BAM_TO_FASTQ_ILLUMINA {
  input:
    tuple val(meta), path(bam)

  output:
  // TODO: Needs to output a path
  // Can it output a pair?
    tuple val(meta), path("${meta.SAMPLE_ID}_*.fastq.gz")

  script:
    fastq_preproc = "fastq_preproc"
    fastq_directory = "fastq_files"

    sorted_bam = "${meta.SAMPLE_ID}.sorted.bam"

    fastq_1_path = "${meta.SAMPLE_ID}_1.fastq.gz"
    fastq_2_path = "${meta.SAMPLE_ID}_2.fastq.gz"
    // reads_map = ["seq_file_1": file(fastq_1_path), "seq_file_2": file(fastq_2_path)]

    """
    # sort by coordinates
    samtools sort -o ${sorted_bam} ${bam}

    # Starting from a coordinate sorted file, output paired reads to separate files, discarding singletons, supplementary and secondary reads. The resulting files can be used with, for example, the bwa aligner.
    # see: http://www.htslib.org/doc/samtools-fasta.html
    samtools collate -u \
      -O ${sorted_bam} \
    | samtools \
      fastq \
      -1 ${fastq_1_path} \
      -2  ${fastq_2_path} \
      -0 /dev/null \
      -s /dev/null \
      -n
    """
}