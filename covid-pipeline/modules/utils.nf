import java.nio.file.Paths


/* return the name of the last dir (e.g. /a/b/c/d => d) */
def getDirName (myDir) {
    pathDir = Paths.get(myDir)
    dirName = pathDir.getFileName().toString()
    return dirName
}

/* make a glob for retrieving illumina fastq files */
def makeFastqSearchPath (illuminaPrefixes, illuminaSuffixes, fastq_exts) {
    def fastqSearchPath = []

    for (suff in illuminaSuffixes) {
        for(ext in fastq_exts) {
            if ( illuminaPrefixes ) {
                for (prefix in illuminaPrefixes) {
                    dirNameGlob = params.directory.replaceAll(/\/+$/, "") + '/' + prefix + suff + ext
                    fastqSearchPath.add(dirNameGlob)
                }
            } else {
                dirNameGlob = params.directory.replaceAll(/\/+$/, "") + '/' + suff + ext
                fastqSearchPath.add(dirNameGlob)
            }
        }
    }

    return fastqSearchPath
}

/* make a glob for retrieving nanopore-medaka fastq files */
def makeNanoporeSearchPath ( ) {
    // note: we need this specific extension
    filePathGlob = params.directory.replaceAll(/\/+$/, "") + '/**' + '/*.fastq'
    return filePathGlob
}

/* make a glob for retrieving bam files for the ncov-illumina workflow */
def makeBamSearchPath ( ) {
    filePathGlob = params.directory.replaceAll(/\/+$/, "") + '/*.bam'
    return filePathGlob
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
    path notification_file, emit: ch_notification_file

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

// Process which acts as a workaround to append boolean value, whether the sample exists in another channel
process append_metadata_match_to_sample_file{
  input:
    tuple val(sample_name), file(sample_file)
    val samples_with_meta

  output:
    tuple val(sample_name), file(sample_file), val(has_meta)

  script:
    has_meta = samples_with_meta.contains(sample_name)

  """
  """
}

// Process which acts as a workaround to append boolean value, whether the sample exists in another channel
process append_metadata_match_to_sample_file_pair{
  input:
    tuple val(sample_name), file(sample_file1), file(sample_file2)
    val samples_with_meta

  output:
    tuple val(sample_name), file(sample_file1), file(sample_file2), val(has_meta)

  script:
    has_meta = samples_with_meta.contains(sample_name)

  """
  """
}

process store_mismatching_files{
    publishDir COVID_PIPELINE_MISSING_METADATA_PATH, mode: 'copy', overwrite: true

    input:
      file(sample_file)

    output:
      file(sample_file)

    script:

    """
    """
}
