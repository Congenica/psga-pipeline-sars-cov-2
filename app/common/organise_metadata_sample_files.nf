include { check_metadata } from './check_metadata.nf'

if ( params.sequencing_technology == "illumina" ) {
    include { stage_sample_file as stage_sample_bam } from './fetch_sample_files.nf'
    include { stage_sample_file_pair as stage_sample_fastq } from './fetch_sample_files.nf'
    include { bam_to_fastq_illumina as bam_to_fastq } from './utils.nf'
} else if ( params.sequencing_technology == "ont" ) {
    include { stage_sample_file as stage_sample_bam } from './fetch_sample_files.nf'
    include { stage_sample_file as stage_sample_fastq } from './fetch_sample_files.nf'
    include { bam_to_fastq_ont as bam_to_fastq } from './utils.nf'
} else if ( params.sequencing_technology == "unknown" ) {
    include { stage_sample_file as stage_sample_fasta } from './fetch_sample_files.nf'
} else {
    throw new Exception("Error: '--sequencing_technology' can only be 'illumina', 'ont' or 'unknown'")
}

/*
 * process the samples.csv and organise sample files in channels.
 */
workflow organise_metadata_sample_files {
    main:
        return
        check_metadata("${params.configPath}samples.csv")
        ch_metadata = check_metadata.out.ch_metadata

        // organise sample input files by file type
        ch_metadata
            .splitCsv(header: true, sep: ',')
            .branch {
                fastq_pair: it.seq_file_2 =~ /\.(fq|fastq?)(?:\.gz)?$/ && it.seq_file_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fastq: it.seq_file_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fasta: it.seq_file_1 =~ /\.(fa|fasta?)(?:\.gz)?$/
                bam: it.seq_file_1 =~ /\.bam$/
            }
            .set { ch_metadata_samples }

        if ( params.sequencing_technology == "unknown" ) {

            // files are FASTA
            ch_sample_files = stage_sample_fasta(ch_metadata_samples.fasta)

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
            ch_sample_files = ch_fastq_files.mix(ch_bam_to_fastq)

        }

    emit:
        ch_metadata
        ch_sample_files
}