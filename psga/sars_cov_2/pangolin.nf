/*
 * Run: pangolin pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process pangolin_pipeline {
  tag "${task.index} - ${reheadered_fasta}"

  input:
    path reheadered_fasta
    /* We parametrise the data directory so that we can change the version
     * of `pangolin-data` without having to release a new version of the
     * pipeline.
     */ 
    path pangolin_data_dir

  output:
    path "${pangolin_out_directory}/${output_filename}", emit: ch_pangolin_lineage_csv

  script:
    sample_name = reheadered_fasta.getSimpleName()
    pangolin_out_directory = "pangolin_output"
    output_filename = "${sample_name}_lineage_report.csv"
  """
  pangolin ${reheadered_fasta} --outdir ${pangolin_out_directory} --outfile ${output_filename} --datadir ${pangolin_data_dir}
  """
}
