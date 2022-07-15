include { pipeline_start } from './pipeline_lifespan.nf'
include { check_metadata } from './check_metadata.nf'
include { fastqc } from './common/fastqc.nf'

if ( params.ncov_workflow == "illumina_artic" ) {
    include { contamination_removal_illumina as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic } from './ncov2019_artic.nf'
    if ( params.filetype == "fastq" ) {
        include { select_sample_file_pair as get_sample_files } from './common/fetch_sample_files.nf'
    } else if ( params.filetype == "bam" ) {
        include { bam_to_fastq } from './common/utils.nf'
        include { select_sample_file as get_sample_files } from './common/fetch_sample_files.nf'
    } else {
        throw new Exception("Error: '--filetype' can only be 'fastq' or 'bam' for 'illumina_artic' workflow")
    }
} else if ( params.ncov_workflow == "medaka_artic" ) {
    include { contamination_removal_ont as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic } from './ncov2019_artic.nf'
    if ( params.filetype == "fastq" ) {
        include { select_sample_file as get_sample_files } from './common/fetch_sample_files.nf'
    } else {
        throw new Exception("Error: '--filetype' can only be 'fastq' for 'medaka_artic' workflow")
    }
} else if ( params.ncov_workflow == "no_ncov" ) {
    if ( params.filetype == "fasta" ) {
        include { select_sample_file as get_sample_files } from './common/fetch_sample_files.nf'
    } else {
        throw new Exception("Error: '--filetype' can only be 'fasta' for 'no_ncov' workflow")
    }
} else {
    throw new Exception("Error: '--ncov_workflow' can only be 'illumina_artic', 'medaka_artic' or 'no_ncov'")
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
            params.ncov_workflow,
            params.filetype,
            params.scheme_repo_url,
            params.scheme_dir,
            params.scheme,
            params.scheme_version
        )

        check_metadata(
            params.metadata,
            params.run,
            params.filetype,
            params.ncov_workflow
        )

        ch_metadata = check_metadata.out.ch_metadata

        if ( params.filetype == "fasta" ) {

            ch_fasta_files = get_sample_files(
                ch_metadata,
                ".fasta"
            )

            // mock ncov
            ch_ncov_qc_csv = Channel.fromPath('/mock_file')

        } else {

            ch_input_files_fastq = Channel.empty()

            if ( params.ncov_workflow == "illumina_artic") {

                if ( params.filetype == "fastq" ) {
                    ch_input_files_fastq = get_sample_files(
                        ch_metadata,
                        ".fastq.gz"
                    )
                } else if ( params.filetype == "bam" ) {
                    ch_input_files_bam = get_sample_files(
                        ch_metadata,
                        ".bam"
                    )
                    ch_input_files_fastq = bam_to_fastq(ch_input_files_bam)
                } else {
                    log.error """\
                        ERROR: illumina workflow supports only fastq or bam file types.
                        Aborting!
                    """
                    System.exit(1)
                }

            } else if ( params.ncov_workflow == "medaka_artic" ) {

                if ( params.filetype == "fastq" ) {
                    ch_input_files_fastq = get_sample_files(
                        ch_metadata,
                        ".fastq"
                    )
                } else {
                    log.error """\
                        ERROR: nanopore / medaka workflow supports only fastq file type.
                        Aborting!
                    """
                    System.exit(1)
                }

            } else {
                log.error """\
                    ERROR: the only supported workflows are illumina_artic and medaka_artic.
                    Aborting!
                """
                System.exit(1)
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
