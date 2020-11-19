#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define ncov-artic-nf pipeline parameters
 */
params.ncov_pipeline_dir = "${HOME}/ps-bahrain-covid/ncov2019-artic-nf"
params.ncov_docker_image = "ncov2019_edited:latest"
params.ncov_prefix = "covid_test"
/* 1 fastq file only for testing */
params.ncov_fastq_sample_dir = "${HOME}/Bahrain_COVID_s3_data_lite/sample_data"
params.ncov_results    = "${HOME}/ncov_results"
params.ncov_done = "ncov.done"
params.python_docker_image = "144563655722.dkr.ecr.eu-west-1.amazonaws.com/congenica/dev/covid-pipeline:python_latest"
params.fasta_storage_dir = "${GENOME_FASTA_PATH}"



log.info """\
COVID pipeline    v 0.0.1
=========================
ncov_pipeline_dir     : $params.ncov_pipeline_dir
ncov_docker_image     : $params.ncov_docker_image
ncov_prefix           : $params.ncov_prefix
ncov_fastq_sample_dir : $params.ncov_fastq_sample_dir
ncov_results          : $params.ncov_results
python_docker_image   : $params.python_docker_image
fasta_storage_dir     : $params.fasta_storage_dir
"""

/*
 * Import modules
 */
include {run_ncov_artic_nf} from './modules.nf'
include {reheader_genome_fasta} from './modules.nf'



workflow {
    run_ncov_artic_nf(
        params.ncov_pipeline_dir,
        params.ncov_docker_image,
        params.ncov_prefix,
        params.ncov_fastq_sample_dir,
        params.ncov_done
    )

    // we flatten the resulting fasta, so that pipeline branches off per-fasta to its own separate processes
    run_ncov_artic_nf.out.ch_fasta_ncov_results \
        .flatten() \
        .set { fasta_to_reheader }
        
    reheadered_fasta = reheader_genome_fasta(fasta_to_reheader)
}
