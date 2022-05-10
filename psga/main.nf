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
include { store_notification as store_metadata_notification } from './common/utils.nf'

include { fastqc } from './common/fastqc.nf'
include { store_fastqc_reports } from './common/fastqc.nf'

if ( params.workflow == "illumina_artic" ) {
    if ( params.filetype == "fastq" ) {
        include { select_sample_file_pair as get_sample_files } from './common/fetch_sample_files.nf'
    } else if ( params.filetype == "bam" ) {
        include { bam_to_fastq } from './common/utils.nf'
        include { select_sample_file as get_sample_files } from './common/fetch_sample_files.nf'
    } else {
        throw new Exception("Error: '--filetype' can only be 'fastq' or 'bam'")
    }
} else if ( params.workflow == "medaka_artic" ) {
    include { select_sample_file as get_sample_files } from './common/fetch_sample_files.nf'
} else {
    throw new Exception("Error: '--workflow' can only be 'illumina_artic' or 'medaka_artic'")
}

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

    store_metadata_notification(
        check_metadata.out.ch_current_session_samples_with_metadata_file
    )

    ch_input_files = Channel.empty()

    if ( params.workflow == "illumina_artic" && params.filetype == "fastq" ) {
        ch_input_files = get_sample_files(
            check_metadata.out.ch_metadata,
            ".fastq.gz"
        )
    } else if ( params.workflow == "illumina_artic" && params.filetype == "bam" ) {
        ch_input_files_prep = get_sample_files(
            check_metadata.out.ch_metadata,
            ".bam"
        )
        ch_input_files = bam_to_fastq(ch_input_files_prep)
    } else if ( params.workflow == "medaka_artic" && params.filetype == "fastq" ) {
        ch_input_files = get_sample_files(
            check_metadata.out.ch_metadata,
            ".fastq"
        )
    } else {
        log.error """\
            ERROR: nanopore / medaka workflow can only run with fastq input files.
            Aborting!
        """
        System.exit(1)
    }

    // run fastqc for all sample files
    fastqc(ch_input_files)
    ch_fastqc_submitted = store_fastqc_reports(
        fastqc.out.ch_fastqc_html_report.collect(),
        fastqc.out.ch_fastqc_zip_report.collect(),
    )

    if ( params.pathogen == "sars_cov_2" ) {
        psga_workflow = sars_cov_2(fastqc.out.ch_input_files)
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
        ch_fastqc_submitted,
        psga_workflow.ch_analysis_run_results_submitted
    )
}
