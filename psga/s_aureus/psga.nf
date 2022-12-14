/*
 * Main workflow for the pathogen: Staphylococcus aureus.
 */

include { fastqc } from './common/fastqc.nf'
include { organise_metadata_sample_files } from './common/organise_metadata_sample_files.nf'
include { bactopia_one } from './bactopia_one.nf'
include { collate_results } from './collate_results.nf'


workflow psga {

    main:
        organise_metadata_sample_files()
        ch_metadata = organise_metadata_sample_files.out.ch_metadata
        ch_sample_files = organise_metadata_sample_files.out.ch_sample_files

        bactopia_one(ch_sample_files, ch_metadata)

        collate_results(
            bactopia_one.out.ch_result_files_json.collect(),
            bactopia_one.out.ch_results_csv.collect()
        )

        ch_analysis_run_results_submitted = collate_results.out.global_csv_file

    emit:
        ch_analysis_run_results_submitted
}
