process store_notification {
  publishDir "${PSGA_OUTPUT_PATH}/notifications", mode: 'copy', overwrite: true

  input:
    path notification_file

  output:
    path notification_file, emit: ch_notification_file

  script:
  """
  """
}
