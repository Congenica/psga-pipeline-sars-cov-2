#!/usr/bin/env nextflow

include { PANGOLIN_PIPELINE } from '../modules/pangolin.nf'

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
        .set { ch_pangolin_input }

    ch_pangolin_input.view()

    reference_data_paths = Channel.fromPath("${params.configPath}reference_data.csv")
        .splitCsv(header: true, 'sep': ',')
        .reduce([]) { result, row -> ["${row.NAME}": row.LOCATION] + result }

    // FASTA
    // [
    //     [SAMPLE_ID:4a32dffd-584c-4f2a-8b17-e1e27b630f66],
    //     /app/work/03/edebd4509e80cedd1a06b058492869/4a32dffd-584c-4f2a-8b17-e1e27b630f66_1.fasta
    // ]
    // Reference data
    // [
    //     medaka:reference/pangolin-data/branch/prerelease/1.23.1/00001/f1b0a3799ebf31edae945d8fab2bd844ecc8952c,
    //     pangolin-data:reference/pangolin-data/branch/prerelease/1.23.1/00001/f1b0a3799ebf31edae945d8fab2bd844ecc8952c
    // ]

    PANGOLIN_PIPELINE(ch_pangolin_input, reference_data_paths)

}