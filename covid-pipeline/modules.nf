/*
 * Run: ncov-artic-nf nextflow pipeline
 */

process run_ncov_artic_nf {
  publishDir params.ncov_results, mode: 'copy', overwrite: true

  input:
    path ncov_pipeline_dir
    val ncov_docker_image
    val ncov_prefix
    path ncov_sample_dir
    val ncov_done

  output:
    path "${ncov_out_directory}/*", emit: all_ncov_results
    path "${ncov_out_directory}/ncovIllumina_sequenceAnalysis_makeConsensus/*.fa", emit: fasta_ncov_results

  script:
    ncov_out_directory = "ncov_output"

  """

  nextflow run ${ncov_pipeline_dir} -profile docker --illumina --prefix ${ncov_prefix} --directory ${ncov_sample_dir} -with-docker ${ncov_docker_image} --outdir ${ncov_out_directory}

  touch ${ncov_done}
  """
}

process reheader_genome_fasta {
  publishDir params.fasta_storage_dir, mode: 'copy', overwrite: true

  container params.python_docker_image

  input:
    path ncov_fasta
    path script_file

  output:
    path "*.fasta"

  script:
    output_directory = "./"

  """
  python /app/scripts/reheader_fasta.py ${ncov_fasta} ${output_directory}
  """
}
