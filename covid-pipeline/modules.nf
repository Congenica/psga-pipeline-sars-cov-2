/*
 * Run: ncov2019-artic-nf nextflow pipeline
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 */
process ncov2019_artic_nf_pipeline {
  publishDir COVID_PIPELINE_WORKDIR, mode: 'copy', overwrite: true

  input:
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
  nextflow run ${COVID_PIPELINE_ROOTDIR}/ncov2019-artic-nf -profile docker --illumina --prefix ${ncov_prefix} --directory ${COVID_PIPELINE_FASTQ_PATH} -with-docker ${ncov_docker_image} --outdir ${ncov_out_directory}
  """
}

process load_ncov_assembly_qc_to_db {
  input:
    tuple val(sample_name), val(pct_N_bases), val(pct_covered_bases), val(longest_no_N_run), val(num_aligned_reads), val(qc_pass)
    file qc_plot_files

  output:
    path ncov_qc_submit_done

  script:
    matching_depth_file_name = "${sample_name}.depth.png"   // unable to map effectively output channel in main, and first and filter functions with regex don't work as well within the process
    ncov_qc_submit_done = "${sample_name}.ncov_qc_submit.done"

  """
  python /app/scripts/submit_sample_qc.py \
      --sample_name ${sample_name} \
      --pct_n_bases ${pct_N_bases} \
      --pct_covered_bases ${pct_covered_bases} \
      --longest_no_n_run ${longest_no_N_run} \
      --num_aligned_reads ${num_aligned_reads} \
      --qc_pass ${qc_pass} \
      --qc_plot ${matching_depth_file_name} \
      --pipeline_version ${workflow.manifest.version}
  
  touch ${ncov_qc_submit_done}
  """
}


/*
 * Reheader the genome fasta files generated by ncov2019-artic-nf
 */
process reheader_genome_fasta {
  publishDir COVID_PIPELINE_FASTA_PATH, mode: 'copy', overwrite: true

  input:
    path ncov_fasta

  output:
    path "*.fasta"

  script:
    files_dir = "./"

  """
  python /app/scripts/reheader_fasta.py ${files_dir}
  """
}

/*
 * Run: pangolin snakemake pipeline
 * see: https://github.com/cov-lineages/pangolin
 */
process pangolin_pipeline {
  publishDir COVID_PIPELINE_WORKDIR, mode: 'copy', overwrite: true

  input:
    path reheadered_fasta

  output:
    path "${pangolin_out_directory}/*", emit: ch_pangolin_results

  script:
    pangolin_out_directory = "pangolin_output"
    output_filename = "${reheadered_fasta.getSimpleName()}_lineage_report.csv"
  """
  pangolin ${reheadered_fasta} --outdir ${pangolin_out_directory} --outfile ${output_filename}
  """
}
