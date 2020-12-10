// Imports from 3rd party pipelines
include { makeFastqSearchPath } from '../ncov2019-artic-nf/modules/util.nf'


workflow filter_fastq_matching_with_metadata{
    take: 
        ch_samples_with_metadata_loaded
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
        // ch_fastq_file_pairs: [[sample_name, fastq1, fastq2], ...]
        // Matching our extracted sample names with samples loaded from I-SEHA
        ch_fastq_file_pair_matches = append_metadata_match_to_sample_file_pair(
            ch_fastq_file_pairs, 
            ch_samples_with_metadata_loaded.collect()
        )

        ch_fastq_file_pair_matches
            .branch { it ->
                matched: it[3]
                    return it
                mismatched: true
                    return it
            }
            .set{ ch_fastq_file_pairs_grouped_by_metadata_search }

        ch_fastq_file_pairs_grouped_by_metadata_search.mismatched.map {
            it -> println("WARNING: Sample ${it[0]} was not found in database with metadata. Files ${it[1]} and ${it[2]} will not be processed")
        }

        ch_fastq_file_pairs_grouped_by_metadata_search.matched
             .map { it -> [ it[1], it[2] ] }
             .flatten()
             .set{ ch_fasta_matching_metadata }

        ch_fastq_file_pairs_grouped_by_metadata_search.mismatched
             .map { it -> [ it[1], it[2] ] }
             .flatten()
             .set{ ch_fasta_mismatching_metadata }

        store_fastas_not_matched(ch_fasta_mismatching_metadata)

        ch_fasta_matching_metadata.ifEmpty {
            println("ERROR: No illumina fastq file pairs found matching samples, provided by I-SEHA sample metdata import. Aborting!")
            System.exit(1)
        }

    emit:
        ch_fasta_matching_metadata
}

// Process which acts as a workaround to append boolean value, whether the sample exists in another channel
process append_metadata_match_to_sample_file_pair{
  input:
    tuple val(sample_name), file(fastq1), file(fastq2)
    val samples_with_meta

  output:
    tuple val(sample_name), file(fastq1), file(fastq2), val(has_meta)

  script:
    has_meta = samples_with_meta.contains(sample_name)

  """
  """
}

process store_fastas_not_matched{
    publishDir COVID_PIPELINE_MISSING_METADATA_PATH, mode: 'copy', overwrite: true

    input:
      file(fastq)
    
    output:
      file(fastq)

    script:

    """
    """
}
