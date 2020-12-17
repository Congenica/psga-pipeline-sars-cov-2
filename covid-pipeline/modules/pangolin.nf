/*
 * Run: pangolin snakemake pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process pangolin_pipeline {
  publishDir "${COVID_PIPELINE_PANGOLIN_PATH}/${workflow.sessionId}", mode: 'copy', overwrite: true

  input:
    path reheadered_fasta

  output:
    tuple val(sample_name), path("${pangolin_out_directory}/${output_filename}"), emit: ch_pangolin_lineage_csv

  script:
    sample_name = reheadered_fasta.getSimpleName()
    pangolin_out_directory = "pangolin_output"
    output_filename = "${sample_name}_lineage_report.csv"
  """
  pangolin ${reheadered_fasta} --outdir ${pangolin_out_directory} --outfile ${output_filename}
  """
}

/*
 * Load pangolin data to the database
 */
process load_pangolin_data_to_db {
  input:
    tuple val(sample_name), file(ch_pangolin_result_csv_file)

  output:
    path ch_pangolin_load_data_done

  script:
    ch_pangolin_load_data_done = "${sample_name}.load_pangolin_data_to_db.done"

  """
  python /app/scripts/load_pangolin_data_to_db.py \
    --pangolin-lineage-report-file "${ch_pangolin_result_csv_file}" \
    --sample-name "${sample_name}"
  touch ${ch_pangolin_load_data_done}
  """
}
