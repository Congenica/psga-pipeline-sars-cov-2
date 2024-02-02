/*
 * Run: pangolin snakemake pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process PANGOLIN_PIPELINE {

  input:
    tuple val(meta), path(fasta)

  output:
    tuple val(meta), path("${pangolin_out_directory}/${output_filename}"), emit: ch_pangolin_lineage_csv

  script:
    sample_id = meta.SAMPLE_ID
    pangolin_out_directory = "pangolin_output"
    output_filename = "${sample_id}_lineage_report.csv"
  """
  pangolin ${reheadered_fasta} \
    --outdir ${pangolin_out_directory} \
    --outfile ${output_filename}
  """
}
