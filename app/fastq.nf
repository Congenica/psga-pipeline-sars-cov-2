#!/usr/bin/env nextflow

// include { BAM_TO_FASTQ_ONT } from './modules/bam_to_fastq.nf'

// process PREPAREONTFASTQ {

//     input:
//     tuple val(meta), path(fastq)

//     // output:
//     // tuple val sample_id path fasta

//     script:
//     """
//     echo "working now!"
//     """
// }

workflow {

    samples = Channel.fromPath("${params.configPath}samples.csv")
        .splitCsv(header: true, sep: ',')
        .map { row ->
            (readKeys, metaKeys) = row.keySet().split { it =~ /^seq_file_/ }
            meta = row.subMap(metaKeys)
            // Stage all reads files as files
            reads = row.subMap(readKeys)
            .collectEntries { key, value ->
                [key, file(value)]
            }
            // Return metadata and reads
            [meta, reads]
        }
        .branch {
                fastq_pair: it[1].seq_file_2 =~ /\.(fq|fastq?)(?:\.gz)?$/ && it[1].seq_file_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fastq: it[1].seq_file_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fasta: it[1].seq_file_1 =~ /\.(fa|fasta?)(?:\.gz)?$/
                bam: it[1].seq_file_1 =~ /\.bam$/
            }

    samples.fastq.view()

    if (params.sequencing_technology == "illumina") {
        // samples.bam -> convert to fastq

        // PREPAREILLUMINAFASTQ()
            // Contamination removal
            // Fastqc
            // autodetection
            // ...ncov..
    } else if (params.sequencing_technology == "ont") {
        // prepareONTFastq()
        BAM_TO_FASTQ_ONT(samples.bam)

        // samples.fastq.mix(BAM_TO_FASTQ_ONT.out)
        // Contamination removal
            // Fastqc
            // autodetection
            // ...ncov..
    } else if (params.sequencing_technology == "unknown" ) {
        // Create FASTA channel
        // Reheader FASTA
    }

    // pangolin
    // .mix

    // Combine results
    // collect

}