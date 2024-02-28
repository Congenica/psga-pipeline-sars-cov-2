#!/usr/bin/env nextflow

include { NCOV2019_ARTIC_NF_PIPELINE } from '../modules/ncov2019_artic.nf'

workflow {

    Channel.fromPath("${params.configPath}samples.csv")
        .splitCsv(header: true, sep: ',')
        .map { row ->
            (readKeys, metaKeys) = row.keySet().split { it =~ /^seq_file_/ }
            meta = row.subMap(metaKeys)
            reads = row.subMap(readKeys).values()
            if (reads.size() == 1) {
                reads = reads[0]
            }
            [meta, reads]
        }
        .set { ch_ncov_input }

    ch_ncov_input.view()

    // ONT
    // [
    //     [SAMPLE_ID:4a32dffd-584c-4f2a-8b17-e1e27b630f66, PRIMER:ARTIC_V3],
    //     /app/work/03/edebd4509e80cedd1a06b058492869/4a32dffd-584c-4f2a-8b17-e1e27b630f66_1.fastq.gz
    // ]

    // ILLUMINA [
    //     [SAMPLE_ID:15d3a409-d84e-4e7f-a832-5980956f8ee3, PRIMER:ARTIC_V3],
    //     [
    //         /app/local_test/illumina/15d3a409-d84e-4e7f-a832-5980956f8ee3_1.fastq,
    //         /app/local_test/illumina/15d3a409-d84e-4e7f-a832-5980956f8ee3_2.fastq
    //     ]
    // ]

    NCOV2019_ARTIC_NF_PIPELINE(ch_ncov_input)

}