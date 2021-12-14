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
    * COVID_PIPELINE_ROOTDIR               : ${COVID_PIPELINE_ROOTDIR}
    * COVID_PIPELINE_BAM_PATH              : ${COVID_PIPELINE_BAM_PATH}
    * COVID_PIPELINE_OUTPUT_DIR            : ${COVID_PIPELINE_OUTPUT_DIR}

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
    COVID_PIPELINE_ROOTDIR,
    COVID_PIPELINE_BAM_PATH,
    COVID_PIPELINE_OUTPUT_DIR
    ]) {
    throw new Exception("Found unset global environment variables. See '[:]' above. Abort")
}


// Import modules
include { load_iseha_metadata } from './modules/iseha_metadata.nf'

include { filter_bam_matching_with_metadata } from './modules/bam_match.nf'

include { bam_to_fastq } from './modules/ncov2019_artic.nf'
include { ncov2019_artic_nf_pipeline } from './modules/ncov2019_artic.nf'
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

    load_iseha_metadata(
        "${COVID_PIPELINE_BAM_PATH}/" + params.metadata_file_name
    )

    load_iseha_metadata.out.ch_all_samples_with_metadata_file
        .splitText().map { it.trim() }.set { ch_all_samples_with_metadata_loaded }
    load_iseha_metadata.out.ch_current_session_samples_with_metadata_file
        .splitText().map { it.trim() }.set { ch_current_session_samples_with_metadata_loaded }
    load_iseha_metadata.out.ch_all_samples_ncov2019_artic_qc_passed_file
        .splitText().map { it.trim() }.set { ch_qc_passed_samples }
    load_iseha_metadata.out.ch_current_session_updated_samples_file
        .splitText().map { it.trim() }.set { ch_updated_samples }

    ch_bam_matching_metadata = filter_bam_matching_with_metadata(
        ch_all_samples_with_metadata_loaded,
        ch_current_session_samples_with_metadata_loaded,
        ch_qc_passed_samples,
        ch_updated_samples
    )

    ch_bam_to_fastq = bam_to_fastq(ch_bam_matching_metadata)

    ncov2019_artic_nf_pipeline(
        ch_bam_to_fastq.collect(),
        params.ncov_prefix
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
        pangolin_pipeline.out.ch_pangolin_lineage_csv
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
