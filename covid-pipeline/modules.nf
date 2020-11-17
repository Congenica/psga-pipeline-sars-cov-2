/*
 * Run: ncov-artic-nf nextflow pipeline
 */

process run_ncov_artic_nf {
  input:
    path ncov_pipeline_dir
    val ncov_docker_image
    val ncov_prefix
    path ncov_sample_dir
    path ncov_results
    val ncov_done

  output:
    path ncov_done, emit: ncov_done_channel

  """
  nextflow run ${ncov_pipeline_dir} -profile docker --illumina --prefix ${ncov_prefix} --directory ${ncov_sample_dir} -with-docker ${ncov_docker_image} --outdir ${ncov_results}

  touch ${ncov_done}
  """
}

