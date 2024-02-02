process FASTQC {
  publishDir "${params.output_path}/fastqc", mode: 'copy', overwrite: true, pattern: '*_fastqc.zip'

  input:
    tuple val(meta), path(read_paths)

  output:
    path "*_fastqc.zip", emit: ch_fastqc_zip_report

  script:
    """
    for fq in "${read_paths}"; do
      fastqc -q --limits /limits.txt --outdir \$PWD \$fq
    done
    """
}
