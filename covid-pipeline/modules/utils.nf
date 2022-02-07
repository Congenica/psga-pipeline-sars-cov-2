import java.nio.file.Paths


def makeNanoporeSearchPath ( ) {

    def nanoporeSearchPath = []

    // Make a glob to recurse directories
    dirNameGlob = params.directory.replaceAll(/\/+$/, "") + '**'

    fileNameGlob = '*.fastq.gz'

    // Build a path
    searchPath = Paths.get(dirNameGlob, fileNameGlob)

    nanoporeSearchPath.add(searchPath.toString())

    return nanoporeSearchPath
}

def makeBamSearchPath ( ) {

    def bamSearchPath = []

    // Make a glob to recurse directories
    dirNameGlob = params.directory.replaceAll(/\/+$/, "") + '**'

    // Make a glob for filenames
    fileNameGlob = '*.bam'

    // Build a path
    searchPath = Paths.get(dirNameGlob, fileNameGlob)

    bamSearchPath.add(searchPath.toString())

    return bamSearchPath
}

process concat_elements_to_single_string{
    input:
      val string_value_list

    output:
      val concatenated_string

    script:
      concatenated_string = string_value_list.join(" ")

    """
    """
}

// Checks, which values in first set have matches in reference set
process append_match_to_values_list{
  input:
    val value_set
    val value_set_reference

  output:
    tuple val(value_set), val(is_matching)

  script:
    is_matching = value_set_reference.contains(value_set)

  """
  """
}

process store_notification_with_values_list{
  publishDir COVID_PIPELINE_NOTIFICATIONS_PATH, mode: 'copy', overwrite: true

  input:
    val file_name
    val values_list

  output:
    file notification_file

  script:
    notification_file = file_name.replace(":", "-")
    content = values_list.join("\n")

  """
touch "${notification_file}"
cat <<EOT >> "${notification_file}"
${content}
EOT
  """
}
