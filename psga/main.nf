#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2


if (params.print_config) {
    include { printMainConfig } from './common/help.nf'
    include { printPathogenConfig } from "./help.nf"
    printMainConfig()
    printPathogenConfig()
    exit 0
}

if (params.help) {
    include { printMainHelp } from './common/help.nf'
    include { printPathogenHelp } from "./help.nf"
    printMainHelp()
    printPathogenHelp()
    exit 0
}

include { psga } from "./psga.nf"


// Required environment variables
// Add new env variables to common/help.nf
if( "[:]" in [
    DOCKER_IMAGE_URI_PATH,
    AWS_CONNECTION_TIMEOUT,
    AWS_MAX_CONNECTIONS,
    AWS_MAX_PARALLEL_TRANSFERS,
    K8S_NODE,
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
if ( params.sequencing_technology == "" ) {
    throw new Exception("Error: '--sequencing_technology' must be defined")
}
if ( params.kit == "" ) {
    throw new Exception("Error: '--kit' must be defined")
}
if ( params.output_path == "" ) {
    throw new Exception("Error: '--output_path' must be defined")
}

workflow {
    psga_workflow = psga()
}
