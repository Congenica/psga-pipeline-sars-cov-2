/*
 * Convert a BAM file to illumina FASTQ (2 reads)
 */
process bam_to_fastq_illumina {
  tag "${task.index} - ${bam}"
  input:
    file bam

  output:
    path "${fastq_directory}/*", emit: ch_bam_to_fastq

  script:
    sample_name = "${bam.baseName}"
    fastq_preproc = "fastq_preproc"
    fastq_directory = "fastq_files"
    // illumina suffix format
    fastq_suffix_1 = "_1.fastq.gz"
    fastq_suffix_2 = "_2.fastq.gz"

  """
  mkdir -p ${fastq_preproc}
  mkdir -p ${fastq_directory}

  # sort by coordinates
  samtools sort -o ${fastq_preproc}/${sample_name}.sorted.bam ${bam}

  # Starting from a coordinate sorted file, output paired reads to separate files, discarding singletons, supplementary and secondary reads. The resulting files can be used with, for example, the bwa aligner.
  # see: http://www.htslib.org/doc/samtools-fasta.html
  samtools collate -u -O ${fastq_preproc}/${sample_name}.sorted.bam | samtools fastq -1 ${fastq_directory}/${sample_name}${fastq_suffix_1} -2  ${fastq_directory}/${sample_name}${fastq_suffix_2} -0 /dev/null -s /dev/null -n
  """
}

/*
 * Convert a BAM file to ONT FASTQ
 */
process bam_to_fastq_ont {
  tag "${task.index} - ${bam}"
  input:
    file bam

  output:
    path "${fastq_output}", emit: ch_bam_to_fastq

  script:
    fastq_output = "${bam.baseName}.fastq"
  """
  # write singleton reads to (decompressed) fastq file. Do not append /1 and 2/ to the read name
  samtools bam2fq -nO ${bam} > ${fastq_output}
  """
}
