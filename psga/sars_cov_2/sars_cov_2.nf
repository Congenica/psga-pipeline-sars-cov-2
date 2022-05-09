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
        ch_fastqc_done
        ch_input_files
        analysis_run
        scheme_repo_url
        scheme_dir
        scheme
        scheme_version

    main:

        ncov2019_artic(
            ch_fastqc_done,
            ch_input_files,
            analysis_run,
            scheme_repo_url,
            scheme_dir,
            scheme,
            scheme_version
        )

        ch_qc_passed_fasta = reheader(
            ncov2019_artic.out.ch_qc_csv_ncov_result,
            ncov2019_artic.out.ch_fasta_ncov_results
        )

        pangolin(ch_qc_passed_fasta)

        ch_analysis_run_results_submitted = submit_analysis_run_results(
            params.run,
            ncov2019_artic.out.ch_qc_csv_ncov_result.collect(),
            ncov2019_artic.out.ch_fasta_ncov_results.collect(),
            ncov2019_artic.out.ch_sample_depth_ncov_results.collect(),
            pangolin.out.ch_pangolin_lineage_csv.collect()
        )

    emit:
        ch_qc_passed_fasta
        ch_analysis_run_results_submitted
}

