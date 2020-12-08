#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
    ======================
    ${workflow.manifest.name} v ${workflow.manifest.version}
    ======================
    ncov2019-artic-nf config:
    * docker image     : $params.ncov_docker_image
    * prefix           : $params.ncov_prefix

    env vars:
    * COVID_PIPELINE_ROOTDIR    : ${COVID_PIPELINE_ROOTDIR}
    * COVID_PIPELINE_FASTQ_PATH : ${COVID_PIPELINE_FASTQ_PATH}
    * COVID_PIPELINE_WORKDIR    : ${COVID_PIPELINE_WORKDIR}
    * COVID_PIPELINE_FASTA_PATH : ${COVID_PIPELINE_FASTA_PATH}
    * DB_HOST                   : ${DB_HOST}
    * DB_NAME                   : ${DB_NAME}
    ======================
"""

// These do not do anything. However, if user environment is missing of these env variables,
// nextflow will not allow to run the pipeline until these env variables are set
required_variable = DB_USER
required_variable = DB_PASSWORD

// Import modules
include { load_iseha_metadata } from './modules.nf'
include { ncov2019_artic_nf_pipeline } from './modules.nf'
include { load_ncov_assembly_qc_to_db } from './modules.nf'
include { prepare_tsv_for_nextstrain } from './modules.nf'
include { reheader_genome_fasta } from './modules.nf'
include { store_reheadered_fasta_passed } from './modules.nf'
include { store_reheadered_fasta_failed } from './modules.nf'
include { concatenate_fasta } from './modules.nf'
include { pangolin_pipeline } from './modules.nf'
include { load_pangolin_data_to_db } from './modules.nf'
include { generate_report_strain_level_and_global_context } from './modules.nf'
include { generate_report_strain_first_seen } from './modules.nf'
include { prepare_microreact_tsv } from './modules.nf'


workflow {

    load_iseha_metadata(
        "${COVID_PIPELINE_FASTQ_PATH}/" + params.metadata_file_name
    )

    ncov2019_artic_nf_pipeline(
        params.ncov_docker_image,
        params.ncov_prefix
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

    pangolin_pipeline(ch_reheadered_fasta)

    ch_pangolin_sample_submitted = load_pangolin_data_to_db(
        pangolin_pipeline.out.ch_pangolin_lineage_csv
    )

    ch_nextstrain_input_tsv = prepare_tsv_for_nextstrain(
        ch_ncov_qc_sample_submitted.collect(),
        ch_pangolin_sample_submitted.collect()
    )

    ch_microreact_input_tsv = prepare_microreact_tsv(
        ch_ncov_qc_sample_submitted.collect(),
        ch_pangolin_sample_submitted.collect()
    )

    Channel
        .fromPath( COVID_PIPELINE_FASTA_PATH )
        .set{ archived_fasta }
    concatenate_fasta(
        params.root_genome_fasta,
        ch_qc_passed_fasta.collect(),
        archived_fasta
    )

    generate_report_strain_level_and_global_context(
        params.pangolearn_lineage_notes_url,
        params.pangolearn_metadata_url,
        params.pangolearn_dir,
        ch_pangolin_sample_submitted.collect(),
    )

    generate_report_strain_first_seen(
        ch_pangolin_sample_submitted.collect(),
    )

}
