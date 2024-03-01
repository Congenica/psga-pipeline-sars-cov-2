/*
 * Run: pangolin snakemake pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process PANGOLIN_PIPELINE {
  // Uncomment to see results locally
  // publishDir "${params.output_path}", mode: 'copy', overwrite: true, pattern: "${pangolin_out_directory}/${output_filename}"

  input:
    tuple val(meta), path(reheadered_fasta)
    val reference_data_map

  output:
    path("${pangolin_out_directory}/${output_filename}"), emit: ch_pangolin_lineage_csv

  script:
    sample_id = meta.SAMPLE_ID
    pangolin_out_directory = "pangolin_output"
    output_filename = "${sample_id}_lineage_report.csv"
    pangolin_data_dir = "${reference_data_map["pangolin-data"]}/pangolin_data/data"

    println("Using reference data from: ${pangolin_data_dir}")

    """
    pangolin ${reheadered_fasta} \
      --outdir ${pangolin_out_directory} \
      --outfile ${output_filename} \
      --datadir ${pangolin_data_dir}
    """
}
