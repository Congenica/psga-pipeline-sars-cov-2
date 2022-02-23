process bam_to_fastq {
  input:
    file bam

  output:
    path "${fastq_directory}/*", emit: ch_bam_to_fastq_files

  script:
    sample_name = "${bam.baseName}"
    fastq_preproc = "fastq_preproc"
    fastq_directory = "fastq_files"
    // illumina suffix format
    fastq_suffix_1 = "S01_L001_R1_001"
    fastq_suffix_2 = "S01_L001_R2_001"

  """
  mkdir -p ${fastq_preproc}
  mkdir -p ${fastq_directory}

  # sort by coordinates
  samtools sort -o ${fastq_preproc}/${sample_name}.sorted.bam ${bam}

  # Starting from a coordinate sorted file, output paired reads to separate files, discarding singletons, supplementary and secondary reads. The resulting files can be used with, for example, the bwa aligner.
  # see: http://www.htslib.org/doc/samtools-fasta.html
  samtools collate -u -O ${fastq_preproc}/${sample_name}.sorted.bam | samtools fastq -1 ${fastq_directory}/${sample_name}_${fastq_suffix_1}.fq -2  ${fastq_directory}/${sample_name}_${fastq_suffix_2}.fq -0 /dev/null -s /dev/null -n

  bgzip ${fastq_directory}/${sample_name}_${fastq_suffix_1}.fq
  bgzip ${fastq_directory}/${sample_name}_${fastq_suffix_2}.fq
  """
}


/*
 * Run: ncov2019-artic-nf nextflow pipeline (illumina workflow)
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 */
process ncov2019_artic_nf_pipeline_illumina {
  input:
    file fastq_file
    val ncov_prefix
    val scheme_repo_url
    val scheme_dir
    val scheme
    val scheme_version

  output:
    path "${ncov_out_directory}/*", emit: ch_all_ncov_results
    path "${ncov_out_directory}/ncovIllumina_sequenceAnalysis_makeConsensus/*.fa", emit: ch_fasta_ncov_results
    path "${ncov_out_directory}/${ncov_prefix}.qc.csv", emit: ch_qc_csv_ncov_result
    path "${ncov_out_directory}/qc_plots/*.depth.png", emit: ch_sample_depth_ncov_results

  script:
    ncov_out_directory = "ncov_output"

  """
  # note: we inject our configuration into ncov to override parameters
  # note: `pwd` is the workdir for this nextflow process
  nextflow run ${COVID_PIPELINE_ROOT_PATH}/ncov2019-artic-nf \
      --illumina \
      --prefix ${ncov_prefix} \
      --directory `eval pwd` \
      --outdir ${ncov_out_directory} \
      --schemeRepoURL ${scheme_repo_url} \
      --schemeDir ${scheme_dir} \
      --scheme ${scheme} \
      --schemeVersion ${scheme_version} \
      -c ${COVID_PIPELINE_ROOT_PATH}/covid-pipeline/ncov-custom.config \
      -c ${COVID_PIPELINE_ROOT_PATH}/covid-pipeline/ncov-illumina-k8s.config
  """
}


/*
 * Run: ncov2019-artic-nf nextflow pipeline (nanopore/medaka workflow)
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 * Note: This runs as a shell block
 */
process ncov2019_artic_nf_pipeline_medaka {
  input:
    file input_dir
    val ncov_prefix
    val scheme_repo_url
    val scheme_dir
    val scheme
    val scheme_version

  output:
    path "ncov_output/*", emit: ch_all_ncov_results
    path "ncov_output/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/*.consensus.fasta", emit: ch_fasta_ncov_results
    path "ncov_output/${ncov_prefix}.qc.csv", emit: ch_qc_csv_ncov_result
    path "ncov_output/qc_plots/*.depth.png", emit: ch_sample_depth_ncov_results

  shell:
  '''

  ncov_out_directory="ncov_output"
  # convert nextflow variables to Bash so that the same format is used
  ncov_prefix=!{ncov_prefix}
  scheme_repo_url=!{scheme_repo_url}
  scheme_dir=!{scheme_dir}
  scheme=!{scheme}
  scheme_version=!{scheme_version}

  # move fastq file to its specific barcode dir
  # these files are located in the nextflow workdir. We need to regenerate the barcode dir
  for fq in *.fastq; do
      barcode="`echo ${fq} | egrep -o 'barcode[[:digit:]]+' | head -n1`"
      mkdir ${barcode}
      mv ${fq} ${barcode}
  done

  # note: we inject our configuration into ncov to override parameters
  # note: `pwd` is the workdir for this nextflow process
  nextflow run ${COVID_PIPELINE_ROOT_PATH}/ncov2019-artic-nf \
      --medaka \
      --prefix ${ncov_prefix} \
      --basecalled_fastq `eval pwd` \
      --outdir ${ncov_out_directory} \
      --schemeRepoURL ${scheme_repo_url} \
      --schemeDir ${scheme_dir} \
      --scheme ${scheme} \
      --schemeVersion ${scheme_version} \
      -c ${COVID_PIPELINE_ROOT_PATH}/covid-pipeline/ncov-custom.config \
      -c ${COVID_PIPELINE_ROOT_PATH}/covid-pipeline/ncov-nanopore-k8s.config


  # this is a code correction to the nanopore medaka workflow in order to restore the correct sample names
  # in file names and file content
  # it is coded here in order to avoid changes to the pipeline
  for input_path in `find barcode*/*.fastq`; do
      barcode=`dirname ${input_path}`
      sample_name=`basename ${input_path%.*}`

      for file_to_update in `find ${ncov_out_directory} -name "${ncov_prefix}_${barcode}*" -type f \\( -name "*.csv" -o -name "*.txt" -o -name "*.fastq" -o -name "*.fasta" \\) -type f`; do
          sed -i "s/${ncov_prefix}_${barcode}/${sample_name}/g" ${file_to_update}
      done

      sed -i "s/${ncov_prefix}_${barcode}/${sample_name}/g" ${ncov_out_directory}/${ncov_prefix}.qc.csv

      for file_to_rename in `find ${ncov_out_directory} -name "${ncov_prefix}_${barcode}*" -type f`; do
          file_dir=`dirname ${file_to_rename}`
          file_name=`basename ${file_to_rename}`
          cd ${file_dir}
          rename "s/${ncov_prefix}_${barcode}/${sample_name}/" ${file_name}
          cd - > /dev/null
      done
  done
  '''
}


/*
 * Store ncov2019_artic output
 * We publish only a single channel. This way multiple channels won't conflict on publish
 */
process store_ncov2019_artic_nf_output {
  publishDir "${COVID_PIPELINE_NCOV_OUTPUT_PATH}/${workflow.sessionId}", mode: 'copy', overwrite: true
  publishDir "${COVID_PIPELINE_NCOV_OUTPUT_PATH}/latest", mode: 'copy', overwrite: true

  input:
    path ch_all_ncov_results

  output:
    path ch_all_ncov_results

  script:

  """
  """
}

/*
 * Load ncov2019-artic-nf data to the database
 */
process load_ncov_data_to_db {
  input:
    file ch_qc_ncov_result_csv_file
    file ch_qc_plot_files
    val ch_analysis_run_name

  output:
    path ch_ncov_qc_load_done, emit: ch_ncov_qc_load_done

  script:
    directory_with_qc_depth_files = "./"
    ch_ncov_qc_load_done = "load_ncov_assembly_qc_to_db.done"

  """
  python /app/scripts/load_ncov_data_to_db.py \
    --ncov-qc-csv-file "${ch_qc_ncov_result_csv_file}" \
    --ncov-qc-depth-directory "${directory_with_qc_depth_files}" \
    --analysis-run-name "${ch_analysis_run_name}"
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
