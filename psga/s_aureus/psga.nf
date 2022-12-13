/*
 * Main workflow for the pathogen: Staphylococcus aureus.
 */

include { fastqc } from './common/fastqc.nf'
include { organise_metadata_sample_files } from './common/organise_metadata_sample_files.nf'
include { bactopia_one } from './bactopia_one.nf'
include { generate_s_aureus_output } from './generate_s_aureus_output.nf'
include { collate_results } from './collate_results.nf'


workflow psga {

    main:
        organise_metadata_sample_files()
        ch_metadata = organise_metadata_sample_files.out.ch_metadata
        ch_sample_files = organise_metadata_sample_files.out.ch_sample_files

        bactopia_one(ch_sample_files, ch_metadata)

        // collect waits for all samples to finish. generate_s_aureus_output only runs once
        generate_s_aureus_output(
            ch_metadata,
            bactopia_one.out.ch_input_files,
            bactopia_one.out.annotation_summary,
            bactopia_one.out.variants_txt_for_csv_file,
            bactopia_one.out.all_software_versions,
            bactopia_one.out.checkm_results_txt,
            bactopia_one.out.assembly_json,
            bactopia_one.out.antimicrobial_protein_report,
            bactopia_one.out.antimicrobial_gene_report,
            bactopia_one.out.mykrobe_versions,
            bactopia_one.out.mlst_versions,
            bactopia_one.out.genome_assembly_name
        )

        collate_results(
            generate_s_aureus_output.out.ch_output_csv_file.collect(),
            generate_s_aureus_output.out.ch_output_json_file.collect()
        )

        ch_analysis_run_results_submitted = collate_results.out.global_csv_file

    emit:
        ch_analysis_run_results_submitted
}
