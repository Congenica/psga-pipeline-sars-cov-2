#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// List of supported pathogens (edit as required)
supported_pathogens = ["sars_cov_2"]

// include help and workflow for all available pathogens
if (supported_pathogens.contains(params.pathogen)) {
    if (params.print_config){
        include {printMainConfig} from './common/help.nf'
        include { printPathogenConfig } from "./${params.pathogen}/help.nf"
        printMainConfig()
        printPathogenConfig()
        exit 0
    }

    if (params.help){
        include {printMainHelp} from './common/help.nf'
        include { printPathogenHelp } from "./${params.pathogen}/help.nf"
        printMainHelp()
        printPathogenHelp()
        exit 0
    }

    include { psga } from "./${params.pathogen}/psga"
} else {
    throw new Exception("Error: unrecognised configuration for pathogen '${params.pathogen}'")
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
