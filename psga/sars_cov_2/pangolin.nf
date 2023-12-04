/*
 * Run: pangolin pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process pangolin_pipeline {
  tag "${task.index} - ${reheadered_fasta}"

  input:
    path reheadered_fasta

  output:
    path "${pangolin_out_directory}/${output_filename}", emit: ch_pangolin_lineage_csv

  script:
    sample_name = reheadered_fasta.getSimpleName()
    pangolin_out_directory = "pangolin_output"
    output_filename = "${sample_name}_lineage_report.csv"
  """
  pangolin ${reheadered_fasta} --outdir ${pangolin_out_directory} --outfile ${output_filename}
  """
}
