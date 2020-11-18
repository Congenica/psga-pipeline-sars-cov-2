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
    path ncov_done, emit: ncov_done_channel
    path "${ncov_out_directory}/*", emit: ncov_results_channel

  script:
    ncov_out_directory = "ncov_output"

  """

  nextflow run ${ncov_pipeline_dir} -profile docker --illumina --prefix ${ncov_prefix} --directory ${ncov_sample_dir} -with-docker ${ncov_docker_image} --outdir ${ncov_out_directory}

  touch ${ncov_done}
  """
}

