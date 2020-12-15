process generate_report_strain_level_and_global_context {
  publishDir COVID_PIPELINE_REPORTS_PATH, mode: 'copy', overwrite: true

  input:
    val pangolearn_lineage_notes_url
    val pangolearn_metadata_url
    val pangolearn_dir
    path pangolin_to_db_submit_completion_flag

  output:
    path output_filename

  script:
    report_name = "strain_level_and_global_context"
    output_filename = "${report_name}_report.csv"

  """
  python /app/scripts/generate_report.py \
    --report "${report_name}" \
    --lineage-notes-url "${pangolearn_lineage_notes_url}" \
    --metadata-url "${pangolearn_metadata_url}" \
    --pangolearn-dir "${pangolearn_dir}" \
    --output "${output_filename}"
  """
}

process generate_report_strain_first_seen {
  publishDir COVID_PIPELINE_REPORTS_PATH, mode: 'copy', overwrite: true

  input:
    path pangolin_to_db_submit_completion_flag

  output:
    path output_filename

  script:
    report_name = "strain_first_seen"
    output_filename = "${report_name}_report.csv"

  """
  python /app/scripts/generate_report.py \
    --report "${report_name}" \
    --output "${output_filename}"
  """
}