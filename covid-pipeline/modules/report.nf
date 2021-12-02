process generate_report_strain_level_and_global_context {
  publishDir COVID_PIPELINE_REPORTS_PATH, mode: 'copy', overwrite: true

  input:
    val pango_designation_lineage_notes_url
    val pango_designation_metadata_url
    val pango_designation_dir
    path pangolin_to_db_submit_completion_flag

  output:
    path output_filename

  script:
    report_name = "strain_level_and_global_context"
    output_filename = "${report_name}_report.csv"

  """
  python /app/scripts/generate_report.py \
    --report "${report_name}" \
    --lineage-notes-url "${pango_designation_lineage_notes_url}" \
    --metadata-url "${pango_designation_metadata_url}" \
    --pango_designation-dir "${pango_designation_dir}" \
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

process generate_report_strain_prevalence {
  publishDir COVID_PIPELINE_REPORTS_PATH, mode: 'copy', overwrite: true

  input:
    path pangolin_to_db_submit_completion_flag

  output:
    path output_filename

  script:
    report_name = "strain_prevalence"
    output_filename = "${report_name}_report.csv"

  """
  python /app/scripts/generate_report.py \
    --report "${report_name}" \
    --output "${output_filename}"
  """
}

process generate_report_sample_dump {
  publishDir COVID_PIPELINE_REPORTS_PATH, mode: 'copy', overwrite: true

  input:
    path pangolin_to_db_submit_completion_flag

  output:
    path output_filename

  script:
    report_name = "sample_dump"
    output_filename = "${report_name}_report.csv"

  """
  python /app/scripts/generate_report.py \
    --report "${report_name}" \
    --output "${output_filename}"
  """
}
