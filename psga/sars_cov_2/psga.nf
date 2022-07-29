include { pipeline_start } from './pipeline_lifespan.nf'
include { check_metadata } from './check_metadata.nf'
include { fastqc } from './common/fastqc.nf'

if ( params.sequencing_technology == "illumina" ) {
    include { contamination_removal_illumina as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic } from './ncov2019_artic.nf'
    include { stage_sample_file } from './common/fetch_sample_files.nf'
    include { stage_sample_file_pair } from './common/fetch_sample_files.nf'
    include { bam_to_fastq } from './common/utils.nf'
} else if ( params.sequencing_technology == "ont" ) {
    include { contamination_removal_ont as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic } from './ncov2019_artic.nf'
    include { stage_sample_file } from './common/fetch_sample_files.nf'
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
    throw new Exception("Found unset environment variables specific to the sars_cov_2 pathoge. See '[:]' above. Abort")
}


/*
 * Main workflow for the pathogen: SARS-CoV-2.
 * This workflow is based on ncov2019-artic and Pangolin pipelines.
 */
workflow psga {

    main:

        // save the session_id and command
        pipeline_start(
            params.metadata,
            params.run,
            params.sequencing_technology,
            params.scheme_repo_url,
            params.scheme_dir,
            params.scheme,
            params.scheme_version
        )

        check_metadata(
            params.metadata,
            params.run,
            params.sequencing_technology
        )
        ch_metadata = check_metadata.out.ch_metadata

        // split the metadata in "single" (1 single file) and "pair" (2 reads) branches
        ch_metadata
            .splitCsv(header: true, sep: ',')
            .branch {
                single: it.file_2 == ''
                pair: true
            }
            .set { ch_metadata_records }

        if ( params.sequencing_technology == "unknown" ) {

            // files are FASTA
            ch_fasta_files = stage_sample_file(ch_metadata_records.single)

            // mock ncov
            ch_ncov_qc_csv = Channel.fromPath('/mock_file')

        } else {

            if ( params.sequencing_technology == "illumina") {

                // stage bam and transform it into 2 paired fastq files
                ch_bam = stage_sample_file(ch_metadata_records.single)
                ch_fastq_pairs_from_bam = bam_to_fastq(ch_bam)

                // stage fastq file pair
                ch_fastq_pairs = stage_sample_file_pair(ch_metadata_records.pair)

                // unify the content of the two channels into one. Order of samples is irrelevant.
                ch_input_files_fastq = ch_fastq_pairs.mix(ch_fastq_pairs_from_bam)

            } else if ( params.sequencing_technology == "ont" ) {

                // stage fastq file pair
                ch_input_files_fastq = stage_sample_file(ch_metadata_records.single)

            } else {

                throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")

            }

            contamination_removal(
                params.rik_ref_genome_fasta,
                ch_input_files_fastq
            )

            fastqc(contamination_removal.out.ch_output_file)

            ncov2019_artic(
                fastqc.out.ch_input_files,
                params.run,
                params.scheme_repo_url,
                params.scheme_dir,
                params.scheme,
                params.scheme_version
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
