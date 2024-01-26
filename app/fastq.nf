#!/usr/bin/env nextflow

include { BAM_TO_FASTQ_ONT } from './modules/bam_to_fastq.nf'
// include { PROCESS_ONT_FASTQ } from './modules/process_fastq.nf'
include { CONTAMINATION_REMOVAL } from './modules/contamination_removal.nf'
include { FASTQC } from './modules/fastqc.nf'
include { PRIMER_AUTODETECTION } from './modules/primer_autodetection.nf'
include { NCOV2019_ARTIC_NF_PIPELINE_MEDAKA } from './modules/ncov2019_artic.nf'

// Params to be moved to config
// TODO: This is copied in by docker. Let's put it not at /
params.rik_ref_genome_fasta = "/MN908947.3.no_poly_A.fa"

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
                fastq_single: it[1].seq_file_1 =~ /\.(fq|fastq?)(?:\.gz)?$/
                fasta: it[1].seq_file_1 =~ /\.(fa|fasta?)(?:\.gz)?$/
                bam: it[1].seq_file_1 =~ /\.bam$/
            }

    samples.fastq_single.view()

    if (params.sequencing_technology == "illumina") {
        // samples.bam -> convert to fastq

        // PREPAREILLUMINAFASTQ()
            // Contamination removal
            // Fastqc
            // autodetection
            // ...ncov..
    } else if (params.sequencing_technology == "ont") {
        BAM_TO_FASTQ_ONT(samples.bam)
        CONTAMINATION_REMOVAL(
            params.rik_ref_genome_fasta,
            samples.fastq_single.mix(BAM_TO_FASTQ_ONT.out)
        )
        FASTQC(
            CONTAMINATION_REMOVAL.out.ch_cleaned_fastq
        )
        PRIMER_AUTODETECTION(
            CONTAMINATION_REMOVAL.out.ch_cleaned_fastq
        )
        // Add the primer as metadata
        ch_ncov_input = PRIMER_AUTODETECTION.out.ch_primer_detected.map {
            it ->
                it[0]["PRIMER"] = it[2].text
                [it[0], it[1]]
        }
        ch_ncov_input.view()

        NCOV2019_ARTIC_NF_PIPELINE_MEDAKA(ch_ncov_input)

    } else if (params.sequencing_technology == "unknown" ) {
        // Create FASTA channel
        // Reheader FASTA
    }

    // pangolin
    // .mix

    // Combine results
    // collect

}