/* This workflow runs per sample */

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
 * Reheader a genome fasta file
 */
process reheader_fasta {
  tag "${task.index} - ${fasta}"
  input:
    path fasta

  output:
    path "*.fasta"

  script:
  """
  python ${PSGA_ROOT_PATH}/scripts/sars_cov_2/reheader_fasta.py --input-dir . --output-dir .
  """
}

/*
 * Process to store fastas, which were marked in ncov pipeline as QC_PASS=TRUE
 */
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

/*
 * Process to store fastas, which were marked in ncov pipeline as QC_PASS=FALSE
 */
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
 * If sequencing_technology is "unknown", ncov is not executed.
 * Return the reheadered passed fasta
 */
workflow reheader {
    take:
       // channel can be a tuple (qc, fasta) or just a fasta file
       ch_sample_fasta
    main:

        if ( params.sequencing_technology == "unknown" ) {
            // ncov was not executed

            ch_reheadered_fasta = reheader_fasta(standardise_fasta(ch_sample_fasta))

            // extend with the sample name [sample_name, fasta]
            // there is no QC, so we assume they all passed
            ch_reheadered_fasta
                .map{ file -> tuple( file.baseName, file ) }
                .set{ ch_reheadered_fasta_qc_passed }

        } else {
            // ncov was executed

            ch_reheadered_fasta = reheader_fasta(ch_sample_fasta)

            // define whether sample is QC_PASSED or QC_FAILED
            // the ncov qc file contains 1 record only
            ch_sample_fasta
                .flatten()
                .filter { /^.*\.qc\.csv$/ }
                .splitCsv(header:true)
                .branch {
                    qc_passed: it.qc_pass =~ /TRUE/
                        return it.sample_name
                    qc_failed: true
                        return it.sample_name
                }
                .set{ ch_sample_name_by_qc }

            // extend with the sample name [sample_name, fasta]
            ch_reheadered_fasta
                .map{ file -> tuple( file.baseName, file ) }
                .set{ ch_reheadered_fasta_with_sample_name }

            // join the tuples by sample_name so that ncov results are organised by QC
            ch_reheadered_fasta_with_sample_name
                .join(ch_sample_name_by_qc.qc_passed)
                .set{ ch_reheadered_fasta_qc_passed }

            ch_reheadered_fasta_with_sample_name
                .join(ch_sample_name_by_qc.qc_failed)
                .set{ ch_reheadered_fasta_qc_failed }

            store_reheadered_qc_failed_fasta(ch_reheadered_fasta_qc_failed)

        }

        ch_qc_passed_fasta = store_reheadered_qc_passed_fasta(ch_reheadered_fasta_qc_passed)

    emit:
        ch_qc_passed_fasta
}
