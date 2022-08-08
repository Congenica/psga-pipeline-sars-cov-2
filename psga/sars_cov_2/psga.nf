include { pipeline_start } from './pipeline_lifespan.nf'
include { check_metadata } from './check_metadata.nf'
include { fastqc } from './common/fastqc.nf'

if ( params.sequencing_technology == "illumina" ) {
    include { contamination_removal_illumina as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic } from './ncov2019_artic.nf'
    include { stage_sample_file } from './common/fetch_sample_files.nf'
    include { stage_sample_file_pair as stage_fastq_sample_file } from './common/fetch_sample_files.nf'
    include { bam_to_fastq_illumina as bam_to_fastq } from './common/utils.nf'
} else if ( params.sequencing_technology == "ont" ) {
    include { contamination_removal_ont as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic } from './ncov2019_artic.nf'
    include { stage_sample_file } from './common/fetch_sample_files.nf'
    include { stage_sample_file as stage_fastq_sample_file } from './common/fetch_sample_files.nf'
    include { bam_to_fastq_ont as bam_to_fastq } from './common/utils.nf'
} else if ( params.sequencing_technology == "unknown" ) {
    include { stage_sample_file } from './common/fetch_sample_files.nf'
} else {
    throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")
}

include { reheader } from './reheader.nf'
include { pangolin_pipeline as pangolin } from './pangolin.nf'
include { submit_analysis_run_results } from './submit_analysis_run_results.nf'


// Required environment variables
if( "[:]" in [
    SARS_COV_2_PIPELINE_DOCKER_IMAGE_TAG,
    NCOV2019_ARTIC_NF_ILLUMINA_DOCKER_IMAGE_TAG,
    NCOV2019_ARTIC_NF_NANOPORE_DOCKER_IMAGE_TAG,
    PANGOLIN_DOCKER_IMAGE_TAG,
    ]) {
    throw new Exception("Found unset environment variables specific to the sars_cov_2 pathogen. See '[:]' above. Abort")
}


primer_scheme_name = "SARS-CoV-2"
if ( params.kit in ['V1', 'V2', 'V3'] ) {
    /* This assignment is required for the ncov2019-artic pipeline when running an ONT workflow.
     * In detail, ncov/ONT expects the primer scheme files to be named after the scheme name.
     * As the scheme name changed from V3 to V4 from 'nCoV-2019' to 'SARS-CoV-2', it happens that
     * ncov/ONT does not find the primer files if these combinations are used:
     * - "SARS-CoV-2" and "V3" (also "V2" and "V1")
     * - "nCoV-2019" and "V4" (and more recent versions)
     * The illumina workflow of ncov2019-artic is not affected by this issue.
     *
     * With this reassignment, we make sure that primer scheme name and file names match.
     */
    primer_scheme_name = 'nCoV-2019'
}


/*
 * Main workflow for the pathogen: SARS-CoV-2.
 * This workflow is based on ncov2019-artic and Pangolin pipelines.
 */
workflow psga {

    main:

        // save the session_id and command
        pipeline_start()

        check_metadata(params.metadata)
        ch_metadata = check_metadata.out.ch_metadata

        // split the metadata in "single" (1 single file) and "pair" (2 reads) branches
        ch_metadata
            .splitCsv(header: true, sep: ',')
            .branch {
                single: it.SEQ_FILE_2 == ''
                pair: true
            }
            .set { ch_metadata_records }


        // organise ch_metadata_records.single by file type
        ch_metadata_records.single
            .filter { it.SEQ_FILE_1 =~ /\.bam$/ }
            .set { ch_metadata_records_single_bam }
        ch_metadata_records.single
            .filter { it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/ }
            .set { ch_metadata_records_single_fastq }
        ch_metadata_records.single
            .filter { it.SEQ_FILE_1 =~ /\.(fa|fasta?)(?:\.gz)?$/ }
            .set { ch_metadata_records_single_fasta }

        ch_metadata_records_pair_fastq = ch_metadata_records.pair


        if ( params.sequencing_technology == "unknown" ) {

            // files are FASTA
            ch_fasta_files = stage_sample_file(ch_metadata_records_single_fasta)
            // mock ncov
            ch_ncov_qc_csv = Channel.fromPath('/mock_file')

        } else {

            // stage bam files
            ch_bam_files = stage_sample_file(ch_metadata_records_single_bam)
            // transform bam into 2 paired fastq files
            ch_bam_to_fastq = bam_to_fastq(ch_bam_files)

            if ( params.sequencing_technology == "illumina") {
                // stage fastq file pair
                ch_fastq_files = stage_fastq_sample_file(ch_metadata_records_pair_fastq)
            } else if ( params.sequencing_technology == "ont" ) {
                // stage fastq file
                ch_fastq_files = stage_fastq_sample_file(ch_metadata_records_single_fastq)
            } else {
                throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")
            }

            // unify the input fastq files with fastq files converted from bam files into one channel. Order of sample tuples is irrelevant.
            ch_input_files_fastq = ch_fastq_files.mix(ch_bam_to_fastq)


            contamination_removal(
                params.rik_ref_genome_fasta,
                ch_input_files_fastq
            )

            fastqc(contamination_removal.out.ch_output_file)

            ncov2019_artic(
                fastqc.out.ch_input_files,
                primer_scheme_name,
            )
            ch_ncov_qc_csv = ncov2019_artic.out.ch_ncov_qc_csv
            ch_fasta_files = ncov2019_artic.out.ch_ncov_sample_fasta
        }

        ch_qc_passed_fasta = reheader(ch_fasta_files)

        pangolin(ch_qc_passed_fasta)

        submit_analysis_run_results(
            ch_metadata,
            ch_ncov_qc_csv.collect(),
            pangolin.out.ch_pangolin_lineage_csv.collect()
        )

        ch_analysis_run_results_submitted = submit_analysis_run_results.out.ch_analysis_run_results_submitted

    emit:
        ch_qc_passed_fasta
        ch_analysis_run_results_submitted
}
