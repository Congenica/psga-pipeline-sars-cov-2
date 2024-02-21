#!/usr/bin/env nextflow

include { BAM_TO_FASTQ_ONT; BAM_TO_FASTQ_ILLUMINA } from './modules/bam_to_fastq.nf'
include { CONTAMINATION_REMOVAL } from './modules/contamination_removal.nf'
include { FASTQC } from './modules/fastqc.nf'
include { PRIMER_AUTODETECTION } from './modules/primer_autodetection.nf'
include { NCOV2019_ARTIC_NF_PIPELINE } from './modules/ncov2019_artic.nf'
include { REHEADER_FASTA } from './modules/reheader_fasta.nf'
include { PANGOLIN_PIPELINE } from './modules/pangolin.nf'
include { SUBMIT_ANALYSIS_RUN_RESULTS } from './modules/submit_analysis_results.nf'

// Params to be moved to config
// TODO: This is copied in by docker. Let's put it not at /
params.rik_ref_genome_fasta = "/MN908947.3.no_poly_A.fa"

workflow {

    ch_metadata = Channel.fromPath("${params.configPath}samples.csv")

    samples = ch_metadata
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

    // TODO: Handle both in one
    // We don't process BAMs right now so for later
    // BAM_TO_FASTQ_ILLUMINA(samples.bam)
    BAM_TO_FASTQ_ONT(samples.bam)

    // ---- Single locally executed workflow option
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
    // Add the primer to metadata
    ch_ncov_input = PRIMER_AUTODETECTION.out.ch_primer_detected.map {
        it ->
            it[0]["PRIMER"] = it[2].text
            [it[0], it[1]]
    }
    // ---- END
    ch_ncov_input.view()

    NCOV2019_ARTIC_NF_PIPELINE(ch_ncov_input)

    // If statement may not be required
    if (params.sequencing_technology == "unknown" ) {
        REHEADER_FASTA(samples.fasta)
    }

    PANGOLIN_PIPELINE(NCOV2019_ARTIC_NF_PIPELINE.out.ch_ncov_sample_fasta.mix(REHEADER_FASTA.out))

    // // submit results
    SUBMIT_ANALYSIS_RUN_RESULTS(
        ch_metadata, // Original metadata input file
        CONTAMINATION_REMOVAL.out.ch_contamination_removal_csv.collect(),
        PRIMER_AUTODETECTION.out.ch_primer_data.collect(),
        NCOV2019_ARTIC_NF_PIPELINE.out.ch_ncov_qc_csv.collect(),
        PANGOLIN_PIPELINE.out.ch_pangolin_lineage_csv.collect(),
    )

}