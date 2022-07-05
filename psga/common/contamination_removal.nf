/*
 * Run: read-it-and-keep
 * https://github.com/GlobalPathogenAnalysisService/read-it-and-keep
 * This tool keeps the reads that match the provided target genome.
 */
process contamination_removal_ont {
  publishDir "${PSGA_OUTPUT_PATH}/contamination_removal", mode: 'copy', overwrite: true, pattern: '*_removed_reads.txt'

  tag "${task.index} - ${file_1}"

  input:
    val ref_genome_fasta
    path file_1

  output:
    path "out/*", emit: ch_output_file
    path "*_removed_reads.txt"

  shell:
  '''
  ref_genome_fasta=!{ref_genome_fasta}
  file_1=!{file_1}
  sample_id=$( echo ${file_1} | cut -d '.' -f1)

  readItAndKeep --tech ont --ref_fasta ${ref_genome_fasta} --reads1 ${file_1} --outprefix out 2>&1 | tee ${sample_id}_removed_reads.txt

  mkdir -p out
  # ncov ont requires reads to be decompressed
  gunzip out.reads.fastq.gz
  mv out.reads.fastq out/${file_1}
  '''
}

process contamination_removal_illumina {
  publishDir "${PSGA_OUTPUT_PATH}/contamination_removal", mode: 'copy', overwrite: true, pattern: '*_removed_reads.txt'

  tag "${task.index} - [${file_1}, ${file_2}]"

  input:
    val ref_genome_fasta
    tuple path(file_1), path(file_2)

  output:
    path "out/*", emit: ch_output_file
    path "*_removed_reads.txt"

  shell:
  '''
  ref_genome_fasta=!{ref_genome_fasta}
  file_1=!{file_1}
  file_2=!{file_2}
  sample_id=$( echo ${file_1} | cut -d '_' -f1 )

  readItAndKeep --tech illumina --ref_fasta ${ref_genome_fasta} --reads1 ${file_1} --reads2 ${file_2} --outprefix out 2>&1 | tee ${sample_id}_removed_reads.txt

  mkdir -p out
  mv out.reads_1.fastq.gz out/${file_1}
  mv out.reads_2.fastq.gz out/${file_2}
  '''
}
