#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define ncov-artic-nf pipeline parameters
 */
params.ncov_pipeline_dir = "${HOME}/ncov2019-artic-nf"
params.ncov_docker_image = "ncov2019_edited:latest"
params.ncov_prefix = "covid_test"
/* 1 fastq file only for testing */
params.ncov_fastq_sample_dir = "${HOME}/Bahrain_COVID_s3_data_lite/sample_data"
params.ncov_results    = "${HOME}/ncov_results"
params.ncov_done = "ncov.done"




log.info """\
COVID pipeline    v 0.0.1
=========================
ncov_pipeline_dir     : $params.ncov_pipeline_dir
ncov_docker_image     : $params.ncov_docker_image
ncov_prefix           : $params.ncov_prefix
ncov_fastq_sample_dir : $params.ncov_fastq_sample_dir
ncov_results          : $params.ncov_results
"""

/*
 * Import modules
 */
include { run_ncov_artic_nf } from './modules.nf'



workflow {
    ncov_done = run_ncov_artic_nf(
        params.ncov_pipeline_dir,
        params.ncov_docker_image,
        params.ncov_prefix,
        params.ncov_fastq_sample_dir,
        params.ncov_results,
        params.ncov_done
    )

}
