/*
 * Main workflow for the pathogen: Staphylococcus aureus.
 */

include { fastqc } from './common/fastqc.nf'
include { organise_metadata_sample_files } from './common/organise_metadata_sample_files.nf'
include { bactopia_one } from './bactopia_one.nf'
include { generate_s_aureus_output } from './generate_s_aureus_output.nf'


workflow psga {

    main:
        organise_metadata_sample_files()
        ch_metadata = organise_metadata_sample_files.out.ch_metadata
        ch_sample_files = organise_metadata_sample_files.out.ch_sample_files

        bactopia_one(ch_sample_files, ch_metadata)

        // collect waits for all samples to finish. generate_s_aureus_output only runs once
        generate_s_aureus_output(
            ch_metadata,
            bactopia_one.out.ch_input_files.collect(),
            bactopia_one.out.annotation_summary.collect(),
            bactopia_one.out.variants_txt_for_csv_file.collect(),
            bactopia_one.out.all_software_versions.collect(),
            bactopia_one.out.checkm_results_txt.collect(),
            bactopia_one.out.assembly_json.collect(),
            bactopia_one.out.antimicrobial_protein_report.collect(),
            bactopia_one.out.antimicrobial_gene_report.collect(),
            bactopia_one.out.mykrobe_versions.collect(),
            bactopia_one.out.mlst_versions.collect()
        )
        ch_analysis_run_results_submitted = generate_s_aureus_output.out.ch_output_csv_file

    emit:
        ch_analysis_run_results_submitted
}
