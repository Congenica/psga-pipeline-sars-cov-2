#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// Import modules
include {printPipelineConfig} from './modules/help.nf'
include {printHelp} from './modules/help.nf'
if (params.print_config){
    printPipelineConfig()
    exit 0
}
if (params.help){
    printHelp()
    exit 0
}

include { pipeline_start } from './modules/pipeline_lifespan.nf'
include { check_metadata } from './modules/check_metadata.nf'
include { store_notification as store_metadata_notification } from './modules/utils.nf'

include { fastqc } from './modules/fastqc.nf'
include { store_fastqc_reports } from './modules/fastqc.nf'

if ( params.workflow == "illumina_artic" ) {
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic_nf_pipeline } from './modules/ncov2019_artic.nf'
    if ( params.filetype == "fastq" ) {
        include { select_sample_file_pair as get_sample_files } from './modules/fetch_sample_files.nf'
    } else if ( params.filetype == "bam" ) {
        include { bam_to_fastq } from './modules/ncov2019_artic.nf'
        include { select_sample_file as get_sample_files } from './modules/fetch_sample_files.nf'
    } else {
        throw new Exception("Error: '--filetype' can only be 'fastq' or 'bam'")
    }
} else if ( params.workflow == "medaka_artic" ) {
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic_nf_pipeline } from './modules/ncov2019_artic.nf'
    include { select_sample_file as get_sample_files } from './modules/fetch_sample_files.nf'
} else {
    throw new Exception("Error: '--workflow' can only be 'illumina_artic' or 'medaka_artic'")
}

include { store_ncov2019_artic_nf_output } from './modules/ncov2019_artic.nf'
include { merge_ncov_qc_files } from './modules/ncov2019_artic.nf'
include { load_ncov_data_to_db } from './modules/ncov2019_artic.nf'
include { reheader_genome_fasta } from './modules/ncov2019_artic.nf'
include { store_reheadered_fasta_passed } from './modules/ncov2019_artic.nf'
include { store_reheadered_fasta_failed } from './modules/ncov2019_artic.nf'

include { pangolin_pipeline } from './modules/pangolin.nf'
include { merge_pangolin_files } from './modules/pangolin.nf'
include { load_pangolin_data_to_db } from './modules/pangolin.nf'

include { create_genbank_submission_files } from './modules/genbank.nf'
include { submit_genbank_files} from './modules/genbank.nf'
include { mark_samples_as_submitted_to_genbank} from './modules/genbank.nf'
include { store_genbank_submission} from './modules/genbank.nf'

include { pipeline_end } from './modules/pipeline_lifespan.nf'


// Required environment variables
// Add new env variables to modules/help.nf
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

    // NCOV2019-ARTIC
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

    ncov2019_artic_nf_pipeline(
        fastqc.out.ch_fastqc_done,
        ch_input_files,
        params.run,
        params.scheme_repo_url,
        params.scheme_dir,
        params.scheme,
        params.scheme_version
    )

    merge_ncov_qc_files(
        ncov2019_artic_nf_pipeline.out.ch_qc_csv_ncov_result.collect()
    )

    store_ncov2019_artic_nf_output(
        ncov2019_artic_nf_pipeline.out.ch_fasta_ncov_results.collect(),
        ncov2019_artic_nf_pipeline.out.ch_sample_depth_ncov_results.collect(),
        merge_ncov_qc_files.out.ch_ncov_qc_all_samples
    )

    ch_ncov_qc_sample_submitted = load_ncov_data_to_db(
        merge_ncov_qc_files.out.ch_ncov_qc_all_samples,
        ncov2019_artic_nf_pipeline.out.ch_sample_depth_ncov_results.collect(),
        params.run
    )

    ch_reheadered_fasta = reheader_genome_fasta(ncov2019_artic_nf_pipeline.out.ch_fasta_ncov_results)


    // Samples are split to QC_PASSED and QC_FAILED
    merge_ncov_qc_files.out.ch_ncov_qc_all_samples
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

    merge_pangolin_files(
        pangolin_pipeline.out.ch_pangolin_lineage_csv.collect()
    )

    ch_pangolin_sample_submitted = load_pangolin_data_to_db(
        merge_pangolin_files.out.ch_pangolin_all_lineages,
        params.run
    )

    Channel
        .fromPath( params.genbank_submission_template )
        .set{ ch_genbank_submission_template }
    create_genbank_submission_files(
        ch_qc_passed_fasta.collect(),
        ch_genbank_submission_template,
        params.genbank_submission_comment,
        params.genbank_submitter_name,
        params.genbank_submitter_account_namespace,
        params.genbank_submission_id_suffix,
        params.run
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
            submit_genbank_files.out.ch_genbank_submission_id,
            params.run
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

    pipeline_end(
        params.run,
        ch_fastqc_submitted,
        ch_ncov_qc_sample_submitted,
        ch_pangolin_sample_submitted
    )
}
