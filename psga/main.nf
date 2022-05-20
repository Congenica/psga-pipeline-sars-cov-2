#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// Import common
include {printPipelineConfig} from './common/help.nf'
include {printHelp} from './common/help.nf'
if (params.print_config){
    printPipelineConfig()
    exit 0
}
if (params.help){
    printHelp()
    exit 0
}

/* Add supported pathogens here */
if ( params.pathogen == "sars_cov_2" ) {
    include { sars_cov_2 as psga } from './sars_cov_2/sars_cov_2.nf'
} else {
    throw new Exception("Error: unrecognised pathogen configuration")
}
include { genbank_submission } from './common/genbank.nf'
include { pipeline_end } from './common/pipeline_lifespan.nf'


// Required environment variables
// Add new env variables to common/help.nf
if( "[:]" in [
    DB_HOST,
    DB_PORT,
    DB_NAME,
    DB_USER,
    DB_PASSWORD,
    PSGA_ROOT_PATH,
    PSGA_OUTPUT_PATH,
    PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH,
    PSGA_CLEANUP_WORKDIR,
    DOCKER_IMAGE_PREFIX,
    PSGA_DOCKER_IMAGE_TAG,
    K8S_PULL_POLICY,
    K8S_SERVICE_ACCOUNT,
    K8S_QUEUE_SIZE,
    K8S_STORAGE_CLAIM_NAME,
    K8S_STORAGE_MOUNT_PATH,
    K8S_PROCESS_MAX_RETRIES,
    K8S_PROCESS_CPU_LOW,
    K8S_PROCESS_CPU_HIGH,
    K8S_PROCESS_MEMORY_VERY_LOW,
    K8S_PROCESS_MEMORY_LOW,
    K8S_PROCESS_MEMORY_MEDIUM,
    K8S_PROCESS_MEMORY_HIGH,
    K8S_PROCESS_MEMORY_VERY_HIGH,
    NXF_WORK,
    NXF_EXECUTOR,
    NXF_ANSI_LOG,
    NXF_OPTS
    ]) {
    throw new Exception("Found unset global environment variables. See '[:]' above. Abort")
}

if ( params.run == "" ) {
    throw new Exception("Error: '--run' must be defined")
}

if ( params.metadata == "" ) {
    throw new Exception("Error: '--metadata' must be defined")
}

workflow {

    psga_workflow = psga()

    genbank_submission(psga_workflow.ch_qc_passed_fasta.collect())

    pipeline_end(
        params.run,
        psga_workflow.ch_analysis_run_results_submitted
    )
}
