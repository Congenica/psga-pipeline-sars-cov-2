/*
 * Make sure that the extension is *.fa
 */
process standardise_fasta {
  tag "${task.index} - ${fasta}"
  input:
    path fasta

  output:
    path "*.consensus.fa"

  shell:
  '''
  fasta=!{fasta}
  mv ${fasta} ${fasta%.*}.consensus.fa
  '''
}

/*
 * Reheader a fasta file
 */
process reheader_fasta {
  tag "${task.index} - ${fasta}"
  input:
    path fasta

  output:
    path "*.fasta"

  script:
  """
  python ${PSGA_ROOT_PATH}/scripts/common/reheader_fasta.py --input-dir . --output-dir .
  """
}

process store_reheadered_qc_passed_fasta {
  tag "${task.index} - ${reheadered_fasta_file}"
  publishDir "${params.output_path}/reheadered-fasta", mode: 'copy', overwrite: true

  input:
    tuple val(sample_name), path(reheadered_fasta_file)

  output:
    path matching_file, emit: ch_qc_passed_fasta

  script:
    matching_file = "${sample_name}.fasta"

  """
  """
}

process store_reheadered_qc_failed_fasta {
  tag "${task.index} - ${reheadered_fasta_file}"
  publishDir "${params.output_path}/reheadered-fasta-qc-failed", mode: 'copy', overwrite: true

  input:
    tuple val(sample_name), path(reheadered_fasta_file)

  output:
    path matching_file, emit: ch_qc_failed_fasta

  script:
    matching_file = "${sample_name}.fasta"

  """
  """
}

/*
 * Reheader a fasta file and store it.
 * If sequencing_technology is "unknown", input samples are fasta
 * Return the reheadered passed fasta
 */
workflow reheader {
    take:
       ch_samples_fasta
       ch_samples_passing_qc
       ch_samples_failing_qc
    main:

        if ( params.sequencing_technology == "unknown" ) {

            ch_reheadered_fasta = reheader_fasta(standardise_fasta(ch_samples_fasta))

            // extend with the sample name [sample_name, fasta]
            // there is no QC, so we assume they all passed
            ch_reheadered_fasta
                .map{ file -> tuple( file.baseName, file ) }
                .set{ ch_reheadered_fasta_qc_passed }

        } else {

            ch_reheadered_fasta = reheader_fasta(ch_samples_fasta)

            // extend with the sample name [sample_name, fasta]
            ch_reheadered_fasta
                .map{ file -> tuple( file.baseName, file ) }
                .set{ ch_reheadered_fasta_with_sample_name }

            // join the tuples by sample_name
            ch_reheadered_fasta_with_sample_name
                .join(ch_samples_passing_qc)
                .set{ ch_reheadered_fasta_qc_passed }

            ch_reheadered_fasta_with_sample_name
                .join(ch_samples_failing_qc)
                .set{ ch_reheadered_fasta_qc_failed }

            store_reheadered_qc_failed_fasta(ch_reheadered_fasta_qc_failed)

        }

        ch_qc_passed_fasta = store_reheadered_qc_passed_fasta(ch_reheadered_fasta_qc_passed)

    emit:
        ch_qc_passed_fasta
}
