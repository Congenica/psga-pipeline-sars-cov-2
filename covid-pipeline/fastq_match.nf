// Imports from 3rd party pipelines
include { makeFastqSearchPath } from '../ncov2019-artic-nf/modules/util.nf'
include { group_fastq_by_metadata_match } from './modules.nf'


workflow filter_fastq_matching_with_metadata{
    take: 
        ch_load_iseha_metadata_completion_flag
    main:
        // Original usage of grouping illumina fastq file pairs and extracting sample name
        ch_fastq_search_path = makeFastqSearchPath( 
            params.illuminaPrefixes, 
            params.illuminaSuffixes, 
            params.fastq_exts 
        )
        Channel.fromFilePairs( ch_fastq_search_path, flat: true)
            .filter{ !( it[0] =~ /Undetermined/ ) }
            .set{ ch_fastq_file_pairs }
        // ch_fastq_file_pairs: [[sample_name, file1, file2], ...]

        // Matching our extracted sample names with samples loaded from I-SEHA
        ch_fastq_file_pair_matches_with_metadata = group_fastq_by_metadata_match(
            ch_fastq_file_pairs,
            ch_load_iseha_metadata_completion_flag
        )
        ch_fastq_file_pair_matches_with_metadata
            .branch {
                matched: it[3] == true
                    return [ it[0], it[1], it[2] ]
                mismatched: true
                    return [ it[0], it[1], it[2] ]
            }
            .set{ ch_fastq_file_pairs_grouped_by_metadata_search }

        ch_fastq_file_pairs_grouped_by_metadata_search.matched
            .map { it -> [ it[1], it[2] ] }
            .flatten()
            .set{ ch_fasta_matching_metadata }

        ch_fasta_matching_metadata.ifEmpty {
            println("ERROR\nNo illumina fastq file pairs found matching samples, provided by I-SEHA sample metdata import. Aborting!")
            System.exit(1)
        }

    emit:
        ch_fasta_matching_metadata
}