#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
    ======================
    ${workflow.manifest.name} v ${workflow.manifest.version}
    ======================
    Global environment variables:
    * DB_HOST                              : ${DB_HOST}
    * DB_PORT                              : ${DB_PORT}
    * DB_NAME                              : ${DB_NAME}
    * DB_USER                              : ${DB_USER}
    * COVID_PIPELINE_ROOT_PATH             : ${COVID_PIPELINE_ROOT_PATH}
    * COVID_PIPELINE_INPUT_PATH            : ${COVID_PIPELINE_INPUT_PATH}
    * COVID_PIPELINE_OUTPUT_PATH           : ${COVID_PIPELINE_OUTPUT_PATH}
    * DOCKER_IMAGE_PREFIX                  : ${DOCKER_IMAGE_PREFIX}
    * DOCKER_IMAGE_TAG                     : ${DOCKER_IMAGE_TAG}
    * K8S_PULL_POLICY                      : ${K8S_PULL_POLICY}
    * K8S_SERVICE_ACCOUNT                  : ${K8S_SERVICE_ACCOUNT}
    * K8S_QUEUE_SIZE                       : ${K8S_QUEUE_SIZE}
    * K8S_STORAGE_CLAIM_NAME               : ${K8S_STORAGE_CLAIM_NAME}
    * K8S_STORAGE_MOUNT_PATH               : ${K8S_STORAGE_MOUNT_PATH}
    * NXF_WORK                             : ${NXF_WORK}
    * NXF_EXECUTOR                         : ${NXF_EXECUTOR}
    * NXF_ANSI_LOG                         : ${NXF_ANSI_LOG}

    Internal environment variables:
    * COVID_PIPELINE_MISSING_METADATA_PATH : ${COVID_PIPELINE_MISSING_METADATA_PATH}
    * COVID_PIPELINE_NCOV_OUTPUT_PATH      : ${COVID_PIPELINE_NCOV_OUTPUT_PATH}
    * COVID_PIPELINE_QC_PLOTS_PATH         : ${COVID_PIPELINE_QC_PLOTS_PATH}
    * COVID_PIPELINE_FASTA_PATH            : ${COVID_PIPELINE_FASTA_PATH}
    * COVID_PIPELINE_FASTA_PATH_QC_FAILED  : ${COVID_PIPELINE_FASTA_PATH_QC_FAILED}
    * COVID_PIPELINE_PANGOLIN_PATH         : ${COVID_PIPELINE_PANGOLIN_PATH}
    * COVID_PIPELINE_GENBANK_PATH          : ${COVID_PIPELINE_GENBANK_PATH}
    * COVID_PIPELINE_NOTIFICATIONS_PATH    : ${COVID_PIPELINE_NOTIFICATIONS_PATH}

    ======================
    params:
    * ncov2019_artic_workflow               : ${params.ncov2019_artic_workflow}
    * input_type                            : ${params.input_type}
    * genbank_submitter_name                : ${params.genbank_submitter_name}
    * genbank_submitter_account_namespace   : ${params.genbank_submitter_account_namespace}
    * genbank_submission_template           : ${params.genbank_submission_template}
    * genbank_storage_remote_url            : ${params.genbank_storage_remote_url}
    * genbank_storage_remote_username       : ${params.genbank_storage_remote_username}
    * genbank_storage_remote_directory      : ${params.genbank_storage_remote_directory}

    ======================
"""

// Required environment variables
if( "[:]" in [
    DB_HOST,
    DB_PORT,
    DB_NAME,
    DB_USER,
    DB_PASSWORD,
    COVID_PIPELINE_ROOT_PATH,
    COVID_PIPELINE_INPUT_PATH,
    COVID_PIPELINE_OUTPUT_PATH,
    ANALYSIS_RUN_NAME,
    DOCKER_IMAGE_PREFIX,
    DOCKER_IMAGE_TAG,
    K8S_PULL_POLICY,
    K8S_SERVICE_ACCOUNT,
    K8S_QUEUE_SIZE,
    K8S_STORAGE_CLAIM_NAME,
    K8S_STORAGE_MOUNT_PATH,
    NXF_WORK,
    NXF_EXECUTOR,
    NXF_ANSI_LOG
    ]) {
    throw new Exception("Found unset global environment variables. See '[:]' above. Abort")
}


// Import modules
include { load_metadata } from './modules/load_metadata.nf'

if ( params.ncov2019_artic_workflow == "illumina" ) {
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic_nf_pipeline } from './modules/ncov2019_artic.nf'
    if ( params.input_type == "fastq" ) {
        include { filter_fastq_matching_with_metadata as filter_input_files_matching_metadata } from './modules/fastq_match.nf'
    } else if ( params.input_type == "bam" ) {
        include { bam_to_fastq } from './modules/ncov2019_artic.nf'
        include { filter_bam_matching_with_metadata as filter_input_files_matching_metadata } from './modules/bam_match.nf'
    } else {
        throw new Exception("Error: input_type can only be 'fastq' or 'bam'")
    }
} else if ( params.ncov2019_artic_workflow == "medaka" ) {
    include { getDirName } from './modules/utils.nf'
    include { filter_nanopore_matching_with_metadata as filter_input_files_matching_metadata } from './modules/nanopore_match.nf'
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic_nf_pipeline } from './modules/ncov2019_artic.nf'
} else {
    throw new Exception("Error: ncov2019_artic_workflow can only be 'illumina' or 'medaka'")
}

include { store_ncov2019_artic_nf_output } from './modules/ncov2019_artic.nf'
include { load_ncov_assembly_qc_to_db } from './modules/ncov2019_artic.nf'
include { reheader_genome_fasta } from './modules/ncov2019_artic.nf'
include { store_reheadered_fasta_passed } from './modules/ncov2019_artic.nf'
include { store_reheadered_fasta_failed } from './modules/ncov2019_artic.nf'
include { store_ncov_qc_plots } from './modules/ncov2019_artic.nf'

include { pangolin_pipeline } from './modules/pangolin.nf'
include { load_pangolin_data_to_db } from './modules/pangolin.nf'

include { create_genbank_submission_files } from './modules/genbank.nf'
include { submit_genbank_files} from './modules/genbank.nf'
include { mark_samples_as_submitted_to_genbank} from './modules/genbank.nf'
include { store_genbank_submission} from './modules/genbank.nf'

include { pipeline_complete } from './modules/pipeline_complete.nf'

workflow {

    // METADATA
    load_metadata(
        "${COVID_PIPELINE_INPUT_PATH}/" + params.metadata_file_name,
        "${ANALYSIS_RUN_NAME}"
    )
    load_metadata.out.ch_all_samples_with_metadata_file
        .splitText().map { it.trim() }.set { ch_all_samples_with_metadata_loaded }
    load_metadata.out.ch_current_session_samples_with_metadata_file
        .splitText().map { it.trim() }.set { ch_current_session_samples_with_metadata_loaded }
    load_metadata.out.ch_all_samples_ncov2019_artic_qc_passed_file
        .splitText().map { it.trim() }.set { ch_qc_passed_samples }
    load_metadata.out.ch_current_session_updated_samples_file
        .splitText().map { it.trim() }.set { ch_updated_samples }


    // NCOV2019-ARTIC
    ch_input_files = Channel.empty()
    ch_input_files_prep = filter_input_files_matching_metadata(
        ch_all_samples_with_metadata_loaded,
        ch_current_session_samples_with_metadata_loaded,
        ch_qc_passed_samples,
        ch_updated_samples
    )

    ncov_prefix = "covid_test"
    if ( params.ncov2019_artic_workflow == "illumina" && params.input_type == "fastq" ) {
        ch_input_files = ch_input_files_prep
    } else if ( params.ncov2019_artic_workflow == "illumina" && params.input_type == "bam" ) {
        ch_input_files = bam_to_fastq(ch_input_files_prep)
    } else if ( params.ncov2019_artic_workflow == "medaka" && params.input_type == "fastq" ) {
        ch_input_files = ch_input_files_prep
        // for ncov nanopore/medaka workflow this is not arbitrary. It must be the name of the full run coming from the lab.
        // This is the name of the input dir
        // e.g. 20200311_1427_X1_FAK72834_a3787181
        ncov_prefix = getDirName("${COVID_PIPELINE_INPUT_PATH}")
    } else {
        log.error """\
            ERROR: nanopore / medaka workflow can only run with fastq input files.
            Aborting!
        """
        System.exit(1)
    }
    ncov2019_artic_nf_pipeline(
        ch_input_files.collect(),
        ncov_prefix
    )


    // Taking only a single output channel and publishing output in separate process after `ncov2019_artic_nf_pipeline`
    // Using single output channel is required to avoid publish conflicts, when two channels attempt to write same file
    store_ncov2019_artic_nf_output(
        ncov2019_artic_nf_pipeline.out.ch_all_ncov_results.collect()
    )

    store_ncov_qc_plots(
        ncov2019_artic_nf_pipeline.out.ch_sample_depth_ncov_results
    )

    ch_ncov_qc_sample_submitted = load_ncov_assembly_qc_to_db(
        ncov2019_artic_nf_pipeline.out.ch_qc_csv_ncov_result,
        ncov2019_artic_nf_pipeline.out.ch_sample_depth_ncov_results
    )

    // flatten so that pipeline branches off by fasta file
    ncov2019_artic_nf_pipeline.out.ch_fasta_ncov_results \
        .flatten() \
        .set { ch_fasta_to_reheader }
    ch_reheadered_fasta = reheader_genome_fasta(ch_fasta_to_reheader)

    // Samples are split to QC_PASSED and QC_FAILED
    ncov2019_artic_nf_pipeline.out.ch_qc_csv_ncov_result
        .splitCsv(header:true)
        .branch {
            qc_passed: it.qc_pass =~ /TRUE/
                return it.sample_name
            qc_failed: true
                return it.sample_name
        }
        .set{ ch_sample_row_by_qc }
    ch_qc_passed_fasta = store_reheadered_fasta_passed(
        ch_reheadered_fasta.collect(),
        ch_sample_row_by_qc.qc_passed.flatten()
    )
    store_reheadered_fasta_failed(
        ch_reheadered_fasta.collect(),
        ch_sample_row_by_qc.qc_failed.flatten()
    )

    pangolin_pipeline(ch_qc_passed_fasta)

    ch_pangolin_sample_submitted = load_pangolin_data_to_db(
        pangolin_pipeline.out.ch_pangolin_lineage_csv,
        "${ANALYSIS_RUN_NAME}"
    )

    Channel
        .fromPath( COVID_PIPELINE_FASTA_PATH )
        .set{ archived_fasta }
    Channel
        .fromPath( params.genbank_submission_template )
        .set{ ch_genbank_submission_template }
    create_genbank_submission_files(
        ch_qc_passed_fasta.collect(),
        archived_fasta,
        ch_genbank_submission_template,
        params.genbank_submission_comment,
        params.genbank_submitter_name,
        params.genbank_submitter_account_namespace,
        params.genbank_submission_id_suffix
    )

    create_genbank_submission_files.out.ch_samples_txt
        .splitText()
        .map { it.trim() }
        .set{ ch_samples_to_submit_to_genbank }

    ch_samples_to_submit_to_genbank
        .ifEmpty{ "NO_SAMPLES" }
        .set {ch_no_samples_flag }

    store_genbank_submission(
        create_genbank_submission_files.out.ch_genbank_xml,
        create_genbank_submission_files.out.ch_genbank_zip,
        create_genbank_submission_files.out.ch_samples_txt,
        create_genbank_submission_files.out.ch_genbank_submission_id
    )

    if ( params.genbank_storage_remote_url && params.genbank_storage_remote_username && params.genbank_storage_remote_password && params.genbank_storage_remote_directory) {
        submit_genbank_files(
            create_genbank_submission_files.out.ch_genbank_xml,
            create_genbank_submission_files.out.ch_genbank_zip,
            create_genbank_submission_files.out.ch_samples_txt,
            create_genbank_submission_files.out.ch_genbank_submission_id,
            ch_no_samples_flag.collect(),
            params.genbank_storage_remote_url,
            params.genbank_storage_remote_username,
            params.genbank_storage_remote_password,
            params.genbank_storage_remote_directory
        )

        mark_samples_as_submitted_to_genbank(
            submit_genbank_files.out.ch_genbank_sample_names_txt,
            ch_no_samples_flag.collect(),
            submit_genbank_files.out.ch_genbank_submission_id
        )
    }
    else {
        log.warn """Missing GenBank upload credentials. Upload to GenBank Submission Portal will be skipped.
            Please set the following parameters in nextflow.config:
                - genbank_submitter_name
                - genbank_submitter_account_namespace
                - genbank_storage_remote_directory
                - genbank_storage_remote_username
                - genbank_storage_remote_password
        """
    }

    pipeline_complete(
        ch_ncov_qc_sample_submitted,
        ch_pangolin_sample_submitted
    )

}
