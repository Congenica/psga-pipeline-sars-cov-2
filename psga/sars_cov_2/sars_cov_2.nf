if ( params.workflow == "illumina_artic" ) {
    include { ncov2019_artic_nf_pipeline_illumina as ncov2019_artic } from './ncov2019_artic.nf'
} else if ( params.workflow == "medaka_artic" ) {
    include { ncov2019_artic_nf_pipeline_medaka as ncov2019_artic } from './ncov2019_artic.nf'
} else {
    throw new Exception("Error: '--workflow' can only be 'illumina_artic' or 'medaka_artic'")
}

include { reheader } from './reheader.nf'
include { pangolin_pipeline as pangolin } from './pangolin.nf'
include { submit_analysis_run_results } from './submit_analysis_run_results.nf'

/*
 * Run the SARS-CoV-2 workflow.
 * This workflow is based on ncov2019-artic and Pangolin pipelines.
 */
workflow sars_cov_2 {

    take:
        ch_input_files

    main:

        ncov2019_artic(
            ch_input_files,
            params.run,
            params.scheme_repo_url,
            params.scheme_dir,
            params.scheme,
            params.scheme_version
        )

        ch_qc_passed_fasta = reheader(ncov2019_artic.out.ch_ncov_sample_fasta)

        pangolin(ch_qc_passed_fasta)

        // flatten the ncov results to make sure we deal with a flat list.
        // E.g. [qc.csv, [fa1, fa2], [png1, png2,..]] => [qc.csv, fa1, fa2, png1, png2]
        // note: flatten() before collect() to execute one single process
        ch_analysis_run_results_submitted = submit_analysis_run_results(
            ncov2019_artic.out.ch_ncov_sample_all_results.flatten().collect(),
            pangolin.out.ch_pangolin_lineage_csv.collect()
        )

    emit:
        ch_qc_passed_fasta
        ch_analysis_run_results_submitted
}

