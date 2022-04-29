import java.nio.file.Paths


/* return the name of the last dir (e.g. /a/b/c/d => d) */
def getDirName (myDir) {
    pathDir = Paths.get(myDir)
    dirName = pathDir.getFileName().toString()
    return dirName
}

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
