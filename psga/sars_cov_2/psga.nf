include { check_metadata } from './check_metadata.nf'
include { fastqc } from './common/fastqc.nf'
include { primer_autodetection } from './common/primer_autodetection.nf'

if ( params.sequencing_technology == "illumina" ) {
    include { contamination_removal_illumina as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic } from './ncov2019_artic.nf'
    include { stage_sample_file as stage_sample_bam } from './common/fetch_sample_files.nf'
    include { stage_sample_file_pair as stage_sample_fastq } from './common/fetch_sample_files.nf'
    include { bam_to_fastq_illumina as bam_to_fastq } from './common/utils.nf'
} else if ( params.sequencing_technology == "ont" ) {
    include { contamination_removal_ont as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic } from './ncov2019_artic.nf'
    include { stage_sample_file as stage_sample_bam } from './common/fetch_sample_files.nf'
    include { stage_sample_file as stage_sample_fastq } from './common/fetch_sample_files.nf'
    include { bam_to_fastq_ont as bam_to_fastq } from './common/utils.nf'
} else if ( params.sequencing_technology == "unknown" ) {
    include { stage_sample_file as stage_sample_fasta } from './common/fetch_sample_files.nf'
} else {
    throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")
}

include { reheader } from './reheader.nf'
include { pangolin_pipeline as pangolin } from './pangolin.nf'
include { submit_analysis_run_results } from './submit_analysis_run_results.nf'


// Required environment variables
if( "[:]" in [
    SARS_COV_2_PIPELINE_DOCKER_IMAGE,
    ]) {
    throw new Exception("Found unset environment variables specific to the sars_cov_2 pathogen. See '[:]' above. Abort")
}


if ( params.sequencing_technology in ["illumina", "ont"] ) {
    /* Store primers in ncov dockerfile with path:
     * PRIMERS (ARTIC, Midnight-IDT, Midnight-ONT, NEB-VarSkip)
     * Store primers in ncov dockerfile with path:
     * /primer_schemes/${scheme}/SARS-CoV-2/${scheme_version}
     * e.g. --kit ARTIC_V4 will use the primers stored in: /primer_schemes/ARTIC/SARS-CoV-2/V4
     */
    split_kit = params.kit.split('_')
    if ( split_kit.length != 2 && params.kit != "unknown" ) {
        throw new Exception("--kit must be either PRIMERNAME_PRIMERVERSION or unknown for illumina/ont samples")
    }
} else if ( params.sequencing_technology == "unknown" ) {
    if ( params.kit != "none" ) {
        throw new Exception("--kit must be 'none' for FASTA samples")
    }
}

/*
 * Main workflow for the pathogen: SARS-CoV-2.
 * This workflow is based on ncov2019-artic and Pangolin pipelines.
 */
workflow psga {

    main:

        check_metadata(params.metadata)
        ch_metadata = check_metadata.out.ch_metadata

        // organise sample input files by file type
        ch_metadata
            .splitCsv(header: true, sep: ',')
            .branch {
                bam: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.bam$/
                fastq: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fastq_pair: it.SEQ_FILE_2 =~ /\.(fq|fastq?)(?:\.gz)?$/ && it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fasta: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.(fa|fasta?)(?:\.gz)?$/
            }
            .set { ch_metadata_samples }

        if ( params.sequencing_technology == "unknown" ) {

            // files are FASTA
            ch_fasta_files = stage_sample_fasta(ch_metadata_samples.fasta)
            // mock channels
            ch_contamination_removal_csv = Channel.empty()
            ch_primer_autodetection_csv = Channel.empty()
            ch_ncov_qc_csv = Channel.empty()

        } else {

            // stage bam files
            ch_bam_files = stage_sample_bam(ch_metadata_samples.bam)
            // transform bam into 2 paired fastq files
            ch_bam_to_fastq = bam_to_fastq(ch_bam_files)

            if ( params.sequencing_technology == "illumina") {
                // stage input file (2 reads)
                ch_fastq_files = stage_sample_fastq(ch_metadata_samples.fastq_pair)
            } else if ( params.sequencing_technology == "ont" ) {
                // stage input file
                ch_fastq_files = stage_sample_fastq(ch_metadata_samples.fastq)
            } else {
                throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")
            }

            // unify the input fastq files with fastq files converted from bam files into one channel. Order of sample tuples is irrelevant.
            ch_input_files_fastq = ch_fastq_files.mix(ch_bam_to_fastq)


            contamination_removal(
                params.rik_ref_genome_fasta,
                ch_input_files_fastq
            )
            ch_contamination_removal_csv = contamination_removal.out.ch_contamination_removal_csv

            fastqc(contamination_removal.out.ch_output_file)

            primer_autodetection(fastqc.out.ch_input_files, "SARS-CoV-2")
            ch_primer_autodetection_csv = primer_autodetection.out.ch_primer_data
            ch_primer_autodetection_files = primer_autodetection.out.ch_files
            /*
             * re-organise ncov input files processed by primer autodetection
             * `it` is:
             * [sample primer txt, fastq] (ont)
             * [sample primer txt, [fastq_1, fastq_2]] (illumina)
             * the code below, returns a tuple of files (e.g. [_,_] or [_,_,_]
             */
            ch_primer_autodetection_files
                .map {
                    it[1] instanceof List ? tuple(it[0], *it[1]) : it
                }
                .set { ch_ncov_input_files }

            ncov2019_artic(ch_ncov_input_files)
            ch_ncov_qc_csv = ncov2019_artic.out.ch_ncov_qc_csv
            ch_fasta_files = ncov2019_artic.out.ch_ncov_sample_fasta
        }

        ch_qc_passed_fasta = reheader(ch_fasta_files)

        pangolin(ch_qc_passed_fasta)

        submit_analysis_run_results(
            ch_metadata,
            ch_contamination_removal_csv.collect(),
            ch_primer_autodetection_csv.collect(),
            ch_ncov_qc_csv.collect(),
            pangolin.out.ch_pangolin_lineage_csv.collect()
        )

        ch_analysis_run_results_submitted = submit_analysis_run_results.out.ch_analysis_run_results_submitted

    emit:
        ch_qc_passed_fasta
        ch_analysis_run_results_submitted
}
