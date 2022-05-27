/*
 * Convert a BAM file into FASTQ (2 reads)
 */
process bam_to_fastq {
  tag "${task.index} - ${bam}"
  input:
    file bam

  output:
    path "${fastq_directory}/*", emit: ch_bam_to_fastq_files

  script:
    sample_name = "${bam.baseName}"
    fastq_preproc = "fastq_preproc"
    fastq_directory = "fastq_files"
    // illumina suffix format
    fastq_suffix_1 = "_1.fastq"
    fastq_suffix_2 = "_2.fastq"

  """
  mkdir -p ${fastq_preproc}
  mkdir -p ${fastq_directory}

  # sort by coordinates
  samtools sort -o ${fastq_preproc}/${sample_name}.sorted.bam ${bam}

  # Starting from a coordinate sorted file, output paired reads to separate files, discarding singletons, supplementary and secondary reads. The resulting files can be used with, for example, the bwa aligner.
  # see: http://www.htslib.org/doc/samtools-fasta.html
  samtools collate -u -O ${fastq_preproc}/${sample_name}.sorted.bam | samtools fastq -1 ${fastq_directory}/${sample_name}${fastq_suffix_1} -2  ${fastq_directory}/${sample_name}${fastq_suffix_2} -0 /dev/null -s /dev/null -n

  bgzip ${fastq_directory}/${sample_name}${fastq_suffix_1}
  bgzip ${fastq_directory}/${sample_name}${fastq_suffix_2}
  """
}
