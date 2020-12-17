/*
 * Run: ncov2019-artic-nf nextflow pipeline
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 */
process ncov2019_artic_nf_pipeline {
  publishDir "${COVID_PIPELINE_NCOV_OUTPUT_PATH}/${workflow.sessionId}", mode: 'copy', overwrite: true

  input:
    file fastq_file
    val ncov_docker_image
    val ncov_prefix

  output:
    path "${ncov_out_directory}/*", emit: ch_all_ncov_results
    path "${ncov_out_directory}/ncovIllumina_sequenceAnalysis_makeConsensus/*.fa", emit: ch_fasta_ncov_results
    path "${ncov_out_directory}/${ncov_prefix}.qc.csv", emit: ch_qc_csv_ncov_result
    path "${ncov_out_directory}/qc_plots/*.depth.png", emit: ch_sample_depth_ncov_results

  script:
    ncov_out_directory = "ncov_output"

  """
  # note: we inject our configuration into ncov to override parameters
  nextflow run ${COVID_PIPELINE_ROOTDIR}/ncov2019-artic-nf -profile docker --illumina --prefix ${ncov_prefix} --directory `eval pwd` -with-docker ${ncov_docker_image} --outdir ${ncov_out_directory} -c ${COVID_PIPELINE_ROOTDIR}/covid-pipeline/ncov-illumina.config
  """
}

/*
 * Load ncov2019-artic-nf assembly qc data to the database
 */
process load_ncov_assembly_qc_to_db {
  input:
    file ch_qc_ncov_result_csv_file
    file ch_qc_plot_files

  output:
    path ch_ncov_qc_load_done, emit: ch_ncov_qc_load_done

  script:
    directory_with_qc_depth_files = "./"
    ch_ncov_qc_load_done = "load_ncov_assembly_qc_to_db.done"

  """
  python /app/scripts/submit_sample_qc.py \
    --ncov_qc_csv_file "${ch_qc_ncov_result_csv_file}" \
    --ncov_qc_depth_directory "${directory_with_qc_depth_files}" \
    --pipeline_version "${workflow.manifest.version}"
  touch ${ch_ncov_qc_load_done}
  """
}

/*
 * Reheader the genome fasta files generated by ncov2019-artic-nf
 */
process reheader_genome_fasta {
  input:
    path ncov_fasta

  output:
    path "*.fasta"

  script:
    files_dir = "./"
    output_dir = "./"

  """
  python /app/scripts/reheader_fasta.py ${files_dir} ${output_dir}
  """
}

/*
 * Process to store fastas, which were marked in ncov pipeline as QC_PASS=TRUE
 */
process store_reheadered_fasta_passed {
  publishDir COVID_PIPELINE_FASTA_PATH, mode: 'copy', overwrite: true

  input:
    path all_reheadered_fasta_files
    val sample_name

  output:
    path matching_file, emit: ch_qc_passed_fasta

  script:
    matching_file = "${sample_name}.fasta"

  """
  """
}

/*
 * Process to store fastas, which were marked in ncov pipeline as QC_PASS=FALSE
 */
process store_reheadered_fasta_failed {
  publishDir COVID_PIPELINE_FASTA_PATH_QC_FAILED, mode: 'copy', overwrite: true

  input:
    path all_reheadered_fasta_files
    val sample_name

  output:
    path matching_file, emit: ch_qc_failed_fasta

  script:
    matching_file = "${sample_name}.fasta"

  """
  """
}

/*
 * Storing .png images of sample qc plots in archive
 */
process store_ncov_qc_plots  {
  publishDir COVID_PIPELINE_QC_PLOTS_PATH, mode: 'copy', overwrite: true

  input:
    path qc_plot_image

  output:
    path qc_plot_image

  """
  """
}