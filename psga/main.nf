#!/usr/bin/env nextflow

import java.nio.file.Files;
import java.nio.file.Paths;

// Enable DSL 2 syntax
nextflow.enable.dsl = 2


// these exceptions would be better handled if Nextflow supported inheritance for configs and workflows.
if (params.pathogen_dir == "") {
    throw new Exception("Pipeline configuration error. Create a pathogen config file which initialises the parameter `pathogen_dir`.")
}

if ( !Files.isDirectory(Paths.get(PSGA_ROOT_PATH, "psga", params.pathogen_dir))
     || Files.notExists(Paths.get(PSGA_ROOT_PATH, "psga", params.pathogen_dir, "help.nf"))
     || Files.notExists(Paths.get(PSGA_ROOT_PATH, "psga", params.pathogen_dir, "psga.nf")) ) {
    throw new Exception("Pipeline configuration error. Create a directory called ${params.pathogen_dir} including the files 'psga.nf' and 'help.nf'.")
}

if (params.print_config) {
    include { printMainConfig } from './common/help.nf'
    include { printPathogenConfig } from "./${params.pathogen_dir}/help.nf"
    printMainConfig()
    printPathogenConfig()
    exit 0
}

if (params.help) {
    include { printMainHelp } from './common/help.nf'
    include { printPathogenHelp } from "./${params.pathogen_dir}/help.nf"
    printMainHelp()
    printPathogenHelp()
    exit 0
}

include { psga } from "./${params.pathogen_dir}/psga.nf"

include { pipeline_end } from './common/pipeline_lifespan.nf'


// Required environment variables
// Add new env variables to common/help.nf
if( "[:]" in [
    PSGA_ROOT_PATH,
    PSGA_OUTPUT_PATH,
    PSGA_INCOMPLETE_ANALYSIS_RUNS_PATH,
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

    pipeline_end(
        params.run,
        psga_workflow.ch_analysis_run_results_submitted
    )
}
