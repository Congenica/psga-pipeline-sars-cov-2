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

include { pipeline_start } from './common/pipeline_lifespan.nf'
include { check_metadata } from './common/check_metadata.nf'
include { store_notification as store_invalid_samples_metadata_notification } from './common/utils.nf'
include { store_notification as store_valid_samples_metadata_notification } from './common/utils.nf'

/* supported pathogens */
if ( params.pathogen == "sars_cov_2" ) {
    include { sars_cov_2 } from './sars_cov_2/sars_cov_2.nf'
} else {
    throw new Exception("Error: '--pathogen' can only be 'sars_cov_2'")
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
    NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG,
    NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG,
    PANGOLIN_DOCKER_IMAGE_TAG,
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

    // save the session_id and command
    pipeline_start(
        params.metadata,
        params.run,
        params.workflow,
        params.filetype,
        params.scheme_repo_url,
        params.scheme_dir,
        params.scheme,
        params.scheme_version
    )

    // METADATA
    check_metadata(
        params.load_missing_samples,
        params.metadata,
        params.run,
        params.scheme,
        params.scheme_version,
        params.filetype,
        params.workflow
    )

    store_valid_samples_metadata_notification(
        check_metadata.out.ch_samples_with_valid_metadata_file
    )
    store_invalid_samples_metadata_notification(
        check_metadata.out.ch_samples_with_invalid_metadata_file
    )

    if ( params.pathogen == "sars_cov_2" ) {
        psga_workflow = sars_cov_2(check_metadata.out.ch_metadata)
    } else {
        log.error """\
            ERROR: Unsupported pathogen.
            Aborting!
        """
        System.exit(1)
    }

    genbank_submission(psga_workflow.ch_qc_passed_fasta.collect())

    pipeline_end(
        params.run,
        psga_workflow.ch_analysis_run_results_submitted
    )
}
