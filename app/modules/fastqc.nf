process FASTQC {
  publishDir "${params.output_path}/fastqc", mode: 'copy', overwrite: true, pattern: '*_fastqc.zip'

  input:
    tuple val(meta), path(reads_file)

  output:
    path "*_fastqc.zip", emit: ch_fastqc_zip_report

  script:
    """
    fastqc -q --limits /limits.txt --outdir \$PWD ${reads_file}
    """
}
