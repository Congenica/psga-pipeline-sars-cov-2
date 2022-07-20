include { fastqc } from './common/fastqc.nf'
include { select_sample_file_pair as get_sample_files } from './common/fetch_sample_files.nf'
include { generate_dummy_output } from './generate_dummy_output.nf'


/*
 * Main workflow for the pathogen: dummy_pathogen.
 */
workflow psga {

    main:

        ch_metadata = Channel.fromPath(params.metadata)

        ch_input_files_fastq = get_sample_files(
            ch_metadata,
            ".fastq.gz"
        )

        fastqc(ch_input_files_fastq)

        generate_dummy_output(
            fastqc.out.ch_input_files.collect()
        )

        ch_analysis_run_results_submitted = generate_dummy_output.out.ch_merged_csv_file

    emit:
        ch_analysis_run_results_submitted
}
