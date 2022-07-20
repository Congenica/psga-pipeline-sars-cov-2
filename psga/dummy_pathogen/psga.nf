include { check_metadata } from './check_metadata.nf'
include { fastqc } from './common/fastqc.nf'

if ( params.sequencing_technology == "illumina" ) {
    include { stage_sample_file as stage_sample_bam } from './common/fetch_sample_files.nf'
    include { stage_sample_file_pair as stage_sample_fastq } from './common/fetch_sample_files.nf'
} else if ( params.sequencing_technology == "ont" ) {
    include { stage_sample_file as stage_sample_bam } from './common/fetch_sample_files.nf'
    include { stage_sample_file as stage_sample_fastq } from './common/fetch_sample_files.nf'
} else {
    throw new Exception("Error: '--sequencing_technology' can only be 'illumina' or 'ont'")
}

include { generate_dummy_output } from './generate_dummy_output.nf'


/*
 * Main workflow for the pathogen: dummy_pathogen.
 */
workflow psga {

    main:
        check_metadata(params.metadata)
        ch_metadata = check_metadata.out.ch_metadata

        // organise sample input files by file type
        ch_metadata
            .splitCsv(header: true, sep: ',')
            .branch {
                fastq: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fastq_pair: it.SEQ_FILE_2 =~ /\.(fq|fastq?)(?:\.gz)?$/ && it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
            }
            .set { ch_metadata_samples }

        if ( params.sequencing_technology == "illumina") {
            // stage input file (2 reads)
            ch_input_files = stage_sample_fastq(ch_metadata_samples.fastq_pair)
        } else if ( params.sequencing_technology == "ont" ) {
            // stage input file
            ch_input_files = stage_sample_fastq(ch_metadata_samples.fastq)
        } else {
                throw new Exception("Error: '--sequencing_technology' can only be 'illumina' or 'ont'")
        }

        fastqc(ch_input_files)

        generate_dummy_output(
            ch_metadata,
            fastqc.out.ch_input_files.collect()
        )

        ch_analysis_run_results_submitted = generate_dummy_output.out.ch_output_csv_file

    emit:
        ch_analysis_run_results_submitted
}
