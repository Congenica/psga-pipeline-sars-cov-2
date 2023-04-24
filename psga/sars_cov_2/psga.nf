include { organise_metadata_sample_files } from './common/organise_metadata_sample_files.nf'
include { fastqc } from './common/fastqc.nf'
include { primer_autodetection } from './common/primer_autodetection.nf'

if ( params.sequencing_technology == "illumina" ) {
    include { contamination_removal_illumina as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic } from './ncov2019_artic.nf'
} else if ( params.sequencing_technology == "ont" ) {
    include { contamination_removal_ont as contamination_removal } from './common/contamination_removal.nf'
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic } from './ncov2019_artic.nf'
} else if ( params.sequencing_technology == "unknown" ) {
    include { reheader_fasta } from './common/reheader.nf'
} else {
    throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")
}

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

        organise_metadata_sample_files()
        ch_metadata = organise_metadata_sample_files.out.ch_metadata
        ch_sample_files = organise_metadata_sample_files.out.ch_sample_files

        if ( params.sequencing_technology == "unknown" ) {

            // files are FASTA
            ch_fasta_files = ch_sample_files
            ch_contamination_removal_csv = Channel.empty()
            ch_primer_autodetection_csv = Channel.empty()
            ch_ncov_qc_csv = Channel.empty()

            ch_reheadered_fasta = reheader_fasta(ch_fasta_files)

        } else {

            // files are FASTQ
            ch_fastq_files = ch_sample_files

            contamination_removal(
                params.rik_ref_genome_fasta,
                ch_fastq_files
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

            // this workflow calls the script reheadering the output fasta file internally for performance reasons.
            // Therefore, the nextflow process reheader_fasta() is not called explicitly
            ncov2019_artic(ch_ncov_input_files)
            ch_ncov_qc_csv = ncov2019_artic.out.ch_ncov_qc_csv
            ch_reheadered_fasta = ncov2019_artic.out.ch_ncov_sample_fasta

        }

        pangolin(ch_reheadered_fasta)

        submit_analysis_run_results(
            ch_metadata,
            ch_contamination_removal_csv.collect(),
            ch_primer_autodetection_csv.collect(),
            ch_ncov_qc_csv.collect(),
            pangolin.out.ch_pangolin_lineage_csv.collect()
        )

        ch_analysis_run_results_submitted = submit_analysis_run_results.out.ch_analysis_run_results_submitted

    emit:
        ch_reheadered_fasta
        ch_analysis_run_results_submitted
}
