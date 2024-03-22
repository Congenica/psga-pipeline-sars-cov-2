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
    pangolin_data_dir = "${reference_data_map["pangolin-data"]}"

    println("Using reference data from: ${pangolin_data_dir}")

    // use-old-datadir is meant to force pangolin to use older reference data
    // even if newer is available.
    // It doesn't seem to do this, so I've kept the older one installed
    // Another alternative would be to install pangolin via pip
    // so that pangolin_data doesn't get installed as a dep by conda

    """
    pangolin ${reheadered_fasta} \
      --outdir ${pangolin_out_directory} \
      --outfile ${output_filename} \
      --datadir ${pangolin_data_dir} \
      --use-old-datadir


    # Add pangolin_data_version column to the end of the result CSV file
    declare -r pangolin_data_version=\$(pangolin --pangolin-data-version | cut -d ' ' -f 2)
    sed -i \
        -e "1s/\$/,pangolin_data_version/" \
        -e "2,\\\$s/\$/,\$pangolin_data_version/" \
        ${pangolin_out_directory}/${output_filename}
    """
}
