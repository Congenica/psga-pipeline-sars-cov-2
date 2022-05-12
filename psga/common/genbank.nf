/*
 * Generate files, required for upload to GenBank
 */
process create_genbank_submission_files {
  input:
    file reheadered_fasta
    file genbank_submission_template
    val genbank_submission_comment
    val genbank_submitter_name
    val genbank_submitter_account_namespace
    val genbank_submission_id_suffix
    val ch_analysis_run_name

  output:
    path submission_xml, emit: ch_genbank_xml
    path submission_zip, emit: ch_genbank_zip
    path submission_samples_txt, emit: ch_samples_txt
    val unique_genbank_submission_identifier, emit: ch_genbank_submission_id

  script:
    sequence_fasta_directory = "./"
    sequence_data_fsa = "sequence.fsa"
    source_metadata_table_src = "source.src"
    submission_xml = "submission.xml"
    submission_zip = "submission.zip"
    submission_samples_txt = "samples_in_submission.txt"
    unique_genbank_submission_identifier = "${workflow.sessionId}.${genbank_submission_id_suffix}"

  """
  echo "All FASTA files to submit to GenBank:"
  ls -l

  python ${PSGA_ROOT_PATH}/scripts/generate_genbank_files.py \
    --analysis-run-name "${ch_analysis_run_name}" \
    --input-sequence-fasta-directory ${sequence_fasta_directory} \
    --input-submission-template ${genbank_submission_template} \
    --output-sequence-data-fsa ${sequence_data_fsa} \
    --output-source-metadata-table-src ${source_metadata_table_src} \
    --output-submission-xml ${submission_xml} \
    --output-submission-zip ${submission_zip} \
    --output-samples-submitted-file ${submission_samples_txt} \
    --submit-name "${genbank_submission_comment}" \
    --submitter "${genbank_submitter_name}" \
    --spuid-namespace "${genbank_submitter_account_namespace}" \
    --spuid-unique-value "${unique_genbank_submission_identifier}"
  """
}

/*
 *  Submit GenBank files to FTP
 */
process submit_genbank_files {
  input:
    path submission_xml
    path submission_zip
    path sample_names_txt
    val submit_id
    val no_samples_flag
    val remote_url
    val remote_username
    val remote_password
    val remote_directory

  when:
    no_samples_flag != ['NO_SAMPLES']

  output:
    path submission_xml, emit: ch_genbank_xml
    path submission_zip, emit: ch_genbank_zip
    path sample_names_txt, emit: ch_genbank_sample_names_txt
    val submit_id, emit: ch_genbank_submission_id

  script:
    ch_submit_genbank_files_done = "submit_genbank_files.done"

  """
  echo ${no_samples_flag}

  python ${PSGA_ROOT_PATH}/scripts/submit_genbank_files.py \
    --input-xml "${submission_xml}" \
    --input-zip "${submission_zip}" \
    --submit-id "${submit_id}" \
    --url "${remote_url}" \
    --username "${remote_username}" \
    --password "${remote_password}" \
    --directory "${remote_directory}"

  touch ${ch_submit_genbank_files_done}
  """
}

/*
 *  Mark samples submitted with GenBank identifier, so that samples are not submitted in next runs
 */
process mark_samples_as_submitted_to_genbank{
  input:
    file sample_names_txt
    val no_samples_flag
    val submit_id
    val ch_analysis_run_name

  when:
    no_samples_flag != ['NO_SAMPLES']

  output:
    file mark_samples_as_submitted_to_genbank_done

  script:
    mark_samples_as_submitted_to_genbank_done = "mark_samples_as_submitted_to_genbank.done"

  """
  python ${PSGA_ROOT_PATH}/scripts/mark_submitted_genbank_samples.py \
    --analysis-run-name "${ch_analysis_run_name}" \
    --sample-names-txt "${sample_names_txt}" \
    --submit-id "${submit_id}"

  touch ${mark_samples_as_submitted_to_genbank_done}
  """
}

/*
 *  Publish GenBank submission in archive directory.
 */
process store_genbank_submission{
  publishDir "${PSGA_OUTPUT_PATH}/genbank", mode: 'copy', overwrite: true

  input:
    path submission_xml
    path submission_zip
    path sample_names_file
    val submit_id

  output:
    path "${out_directory}/*", emit: ch_genbank_submit_files

  script:
    out_directory = submit_id

  """
  mkdir "${out_directory}"
  mv "${submission_xml}" "${out_directory}/${submission_xml}"
  mv "${submission_zip}" "${out_directory}/${submission_zip}"
  mv "${sample_names_file}" "${out_directory}/${sample_names_file}"
  """
}


workflow genbank_submission {
    take:
        ch_qc_passed_fastas
    main:

        Channel
            .fromPath( params.genbank_submission_template )
            .set{ ch_genbank_submission_template }

        create_genbank_submission_files(
            ch_qc_passed_fastas,
            ch_genbank_submission_template,
            params.genbank_submission_comment,
            params.genbank_submitter_name,
            params.genbank_submitter_account_namespace,
            params.genbank_submission_id_suffix,
            params.run
        )

        create_genbank_submission_files.out.ch_samples_txt
            .splitText()
            .map { it.trim() }
            .set{ ch_samples_to_submit_to_genbank }

        ch_samples_to_submit_to_genbank
            .ifEmpty{ "NO_SAMPLES" }
            .set {ch_no_samples_flag }

        store_genbank_submission(
            create_genbank_submission_files.out.ch_genbank_xml,
            create_genbank_submission_files.out.ch_genbank_zip,
            create_genbank_submission_files.out.ch_samples_txt,
            create_genbank_submission_files.out.ch_genbank_submission_id
        )

        if ( params.genbank_storage_remote_url
          && params.genbank_storage_remote_username
          && params.genbank_storage_remote_password
          && params.genbank_storage_remote_directory) {

            submit_genbank_files(
                create_genbank_submission_files.out.ch_genbank_xml,
                create_genbank_submission_files.out.ch_genbank_zip,
                create_genbank_submission_files.out.ch_samples_txt,
                create_genbank_submission_files.out.ch_genbank_submission_id,
                ch_no_samples_flag.collect(),
                params.genbank_storage_remote_url,
                params.genbank_storage_remote_username,
                params.genbank_storage_remote_password,
                params.genbank_storage_remote_directory
            )

            mark_samples_as_submitted_to_genbank(
                submit_genbank_files.out.ch_genbank_sample_names_txt,
                ch_no_samples_flag.collect(),
                submit_genbank_files.out.ch_genbank_submission_id,
                params.run
            )
        }
        else {
            log.warn """Missing GenBank upload credentials. Upload to GenBank Submission Portal will be skipped.
                Please set the following parameters in nextflow.config:
                    - genbank_submitter_name
                    - genbank_submitter_account_namespace
                    - genbank_storage_remote_directory
                    - genbank_storage_remote_username
                    - genbank_storage_remote_password
            """
        }

}
