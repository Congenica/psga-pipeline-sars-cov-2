process BAM_TO_FASTQ {
  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path(output_path)

  script:

  if( params.sequencing_technology == 'illumina' ) {
    sorted_bam = "${meta.SAMPLE_ID}.sorted.bam"

    fastq_1_path = "${meta.SAMPLE_ID}_1.fastq.gz"
    fastq_2_path = "${meta.SAMPLE_ID}_2.fastq.gz"
    output_path = "${meta.SAMPLE_ID}_*.fastq.gz"

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

  } else if( params.sequencing_technology == 'ont' ) {

    output_path = "${meta.SAMPLE_ID}_1.fastq"
    """
    # write singleton reads to (decompressed) fastq file. Do not append /1 and 2/ to the read name
    samtools bam2fq -nO ${bam} > ${output_path}

    """
  }
}