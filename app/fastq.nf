#!/usr/bin/env nextflow

include { BAM_TO_FASTQ_ONT; BAM_TO_FASTQ_ILLUMINA } from './modules/bam_to_fastq.nf'
// include { PROCESS_ONT_FASTQ } from './modules/process_fastq.nf'
include { CONTAMINATION_REMOVAL } from './modules/contamination_removal.nf'
include { FASTQC } from './modules/fastqc.nf'
include { PRIMER_AUTODETECTION } from './modules/primer_autodetection.nf'
include { NCOV2019_ARTIC_NF_PIPELINE } from './modules/ncov2019_artic.nf'

// Params to be moved to config
// TODO: This is copied in by docker. Let's put it not at /
params.rik_ref_genome_fasta = "/MN908947.3.no_poly_A.fa"

workflow {

    samples = Channel.fromPath("${params.configPath}samples.csv")
        .splitCsv(header: true, sep: ',')
        .map { row ->
            (readKeys, metaKeys) = row.keySet().split { it =~ /^seq_file_/ }
            meta = row.subMap(metaKeys)
            reads = row.subMap(readKeys).values()
            [meta, reads]
        }
        .branch {
            fastq: it[1][0] =~ /\.(fq|fastq?)(?:\.gz)?$/
            fasta: it[1][0] =~ /\.(fa|fasta?)(?:\.gz)?$/
            bam: it[1][0] =~ /\.bam$/
            }

    samples.fastq.view()

    // TODO: Make generic
    // BAM_TO_FASTQ_ILLUMINA(samples.bam)
    BAM_TO_FASTQ_ONT(samples.bam)
    CONTAMINATION_REMOVAL(
        params.rik_ref_genome_fasta,
        samples.fastq.mix(BAM_TO_FASTQ_ONT.out)
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

    NCOV2019_ARTIC_NF_PIPELINE(ch_ncov_input)

    if (params.sequencing_technology == "unknown" ) {
        // Create FASTA channel
        // Reheader FASTA
    }

    // pangolin
    // .mix

    // Combine results
    // collect

}