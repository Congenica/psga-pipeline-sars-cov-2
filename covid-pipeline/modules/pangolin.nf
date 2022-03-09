/*
 * Run: pangolin snakemake pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process pangolin_pipeline {

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

/*
 * Merge pangolin lineages into one single file
 */
process merge_pangolin_files {
  publishDir "${COVID_PIPELINE_PANGOLIN_PATH}/${params.run}", mode: 'copy', overwrite: true
  input:
    file input_dir

  output:
    path "${output_path}", emit: ch_pangolin_all_lineages

  script:
    output_path = "all_lineages_report.csv"

  """
  # extract the header from a lineage report
  a_lineage_report="\$(ls *_lineage_report.csv | head -n 1)"
  # copy over the header only (first line)
  awk 'FNR == 1' "\${a_lineage_report}" > ${output_path}
  # copy over the record only from all files (second line)
  awk 'FNR == 2' *_lineage_report.csv >> ${output_path}
  """
}

/*
 * Load pangolin data to the database
 */
process load_pangolin_data_to_db {
  input:
    file ch_pangolin_all_lineages
    val ch_analysis_run_name

  output:
    path ch_pangolin_load_data_done

  script:
    ch_pangolin_load_data_done = "load_pangolin_data_to_db.done"

  """
  python /app/scripts/load_pangolin_data_to_db.py \
    --pangolin-lineage-report-file "${ch_pangolin_all_lineages}" \
    --analysis-run-name "${ch_analysis_run_name}"
  touch ${ch_pangolin_load_data_done}
  """
}
