include { organise_metadata_sample_files } from './common/organise_metadata_sample_files.nf'
include { fastqc } from './common/fastqc.nf'
include { generate_synthetic_output } from './generate_synthetic_output.nf'


/*
 * Main workflow for the pathogen: dummy_pathogen.
 */
workflow psga {

    main:
        organise_metadata_sample_files()
        ch_metadata = organise_metadata_sample_files.out.ch_metadata
        ch_sample_files = organise_metadata_sample_files.out.ch_sample_files

        fastqc(ch_sample_files)

        generate_synthetic_output(
            ch_metadata,
            fastqc.out.ch_input_files.collect()
        )

        ch_analysis_run_results_submitted = generate_synthetic_output.out.ch_output_csv_file

    emit:
        ch_analysis_run_results_submitted
}
