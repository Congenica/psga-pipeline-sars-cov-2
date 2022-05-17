/* This workflow runs per sample */

/*
 * Reheader the genome fasta files generated by ncov2019-artic-nf
 */
process reheader_fasta {
  tag "${task.index} - ${ncov_fasta}"
  input:
    path ncov_fasta

  output:
    path "*.fasta"

  script:
    files_dir = "./"
    output_dir = "./"

  """
  python ${PSGA_ROOT_PATH}/scripts/reheader_fasta.py ${files_dir} ${output_dir}
  """
}

/*
 * Process to store fastas, which were marked in ncov pipeline as QC_PASS=TRUE
 */
process store_ncov2019_artic_qc_passed_fasta {
  tag "${task.index} - ${reheadered_fasta_file}"
  publishDir "${PSGA_OUTPUT_PATH}/reheadered-fasta", mode: 'copy', overwrite: true

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
process store_ncov2019_artic_qc_failed_fasta {
  tag "${task.index} - ${reheadered_fasta_file}"
  publishDir "${PSGA_OUTPUT_PATH}/reheadered-fasta-qc-failed", mode: 'copy', overwrite: true

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
 * Reheader a fasta file and store it based on ncov QC.
 * Return the reheadered fasta with QC status: PASSED
 */
workflow reheader {
    take:
       // a tuple channel: (qc, fasta)
       ch_ncov_sample_fasta
    main:

        ch_reheadered_fasta = reheader_fasta(ch_ncov_sample_fasta)

        // define whether sample is QC_PASSED or QC_FAILED
        // the ncov qc file contains 1 record only
        ch_ncov_sample_fasta
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

        // extend the ncov sample results channel with the sample name
        // [sample_name, fasta]
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


        ch_qc_passed_fasta = store_ncov2019_artic_qc_passed_fasta(ch_reheadered_fasta_qc_passed)
        store_ncov2019_artic_qc_failed_fasta(ch_reheadered_fasta_qc_failed)

    emit:
        ch_qc_passed_fasta
}