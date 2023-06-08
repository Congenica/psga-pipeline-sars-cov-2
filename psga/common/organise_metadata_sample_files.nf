include { check_metadata } from './check_metadata.nf'


/*
 * process the metadata.csv and organise sample files in channels.
 */
workflow organise_metadata_sample_files {
    main:
        check_metadata(params.metadata)
        ch_metadata = check_metadata.out.ch_metadata

        // organise sample input files by file type
        ch_metadata
            .splitCsv(header: true, sep: ',')
            .branch {
                bam: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.bam$/
                    return file(row.SEQ_FILE_1)
                fastq: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                    return file(row.SEQ_FILE_1)
                fastq_pair: it.SEQ_FILE_2 =~ /\.(fq|fastq?)(?:\.gz)?$/ && it.SEQ_FILE_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                    return tuple(file(row.SEQ_FILE_1), file(row.SEQ_FILE_2))
                fasta: it.SEQ_FILE_2 == '' && it.SEQ_FILE_1 =~ /\.(fa|fasta?)(?:\.gz)?$/
                    return file(row.SEQ_FILE_1)
            }
            .set { ch_metadata_samples }

        if ( params.sequencing_technology == "unknown" ) {

            // files are FASTA
            ch_sample_files = ch_metadata_samples.fasta

        } else {

            // transform bam into 2 paired fastq files
            ch_bam_to_fastq = bam_to_fastq(ch_metadata_samples.bam)

            if ( params.sequencing_technology == "illumina") {
                ch_fastq_files = ch_metadata_samples.fastq_pair
            } else if ( params.sequencing_technology == "ont" ) {
                ch_fastq_files = ch_metadata_samples.fastq
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
