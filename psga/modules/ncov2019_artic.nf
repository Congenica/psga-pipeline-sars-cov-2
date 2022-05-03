process bam_to_fastq {
  tag "${task.index} - ${bam}"
  input:
    file bam

  output:
    path "${fastq_directory}/*", emit: ch_bam_to_fastq_files

  script:
    sample_name = "${bam.baseName}"
    fastq_preproc = "fastq_preproc"
    fastq_directory = "fastq_files"
    // illumina suffix format
    fastq_suffix_1 = "_1.fastq"
    fastq_suffix_2 = "_2.fastq"

  """
  mkdir -p ${fastq_preproc}
  mkdir -p ${fastq_directory}

  # sort by coordinates
  samtools sort -o ${fastq_preproc}/${sample_name}.sorted.bam ${bam}

  # Starting from a coordinate sorted file, output paired reads to separate files, discarding singletons, supplementary and secondary reads. The resulting files can be used with, for example, the bwa aligner.
  # see: http://www.htslib.org/doc/samtools-fasta.html
  samtools collate -u -O ${fastq_preproc}/${sample_name}.sorted.bam | samtools fastq -1 ${fastq_directory}/${sample_name}${fastq_suffix_1} -2  ${fastq_directory}/${sample_name}${fastq_suffix_2} -0 /dev/null -s /dev/null -n

  bgzip ${fastq_directory}/${sample_name}${fastq_suffix_1}
  bgzip ${fastq_directory}/${sample_name}${fastq_suffix_2}
  """
}


/*
 * Run: ncov2019-artic-nf nextflow pipeline (illumina workflow)
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 */
process ncov2019_artic_nf_pipeline_illumina {
  tag "${task.index} - ${fastq_file}"
  input:
    file fastq_file
    val ncov_prefix
    val scheme_repo_url
    val scheme_dir
    val scheme
    val scheme_version

  output:
    path "ncov_output/ncovIllumina_sequenceAnalysis_makeConsensus/*.fa", emit: ch_fasta_ncov_results
    path "ncov_output/*.qc.csv", emit: ch_qc_csv_ncov_result
    path "ncov_output/qc_plots/*.png", emit: ch_sample_depth_ncov_results

  shell:
  '''
  ncov_out_directory="ncov_output"

  # convert nextflow variables to Bash so that the same format is used
  ncov_prefix=!{ncov_prefix}
  scheme_repo_url=!{scheme_repo_url}
  scheme_dir=!{scheme_dir}
  scheme=!{scheme}
  scheme_version=!{scheme_version}

  # note: we inject our configuration into ncov to override parameters
  # note: `pwd` is the workdir for this nextflow process
  nextflow run ${PSGA_ROOT_PATH}/ncov2019-artic-nf \
      --illumina \
      --prefix ${ncov_prefix} \
      --directory `eval pwd` \
      --outdir ${ncov_out_directory} \
      --schemeRepoURL ${scheme_repo_url} \
      --schemeDir ${scheme_dir} \
      --scheme ${scheme} \
      --schemeVersion ${scheme_version} \
      -c ${PSGA_ROOT_PATH}/psga/ncov-custom.config \
      -c ${PSGA_ROOT_PATH}/psga/ncov-illumina-k8s.config

  # extract the sample name from one of the reads and rename the qc csv
  sample_name="$(ls *.gz | head -n 1 | cut -d"_" -f1)"
  mv ${ncov_out_directory}/${ncov_prefix}.qc.csv ${ncov_out_directory}/${sample_name}.qc.csv
  '''
}


/*
 * Run: ncov2019-artic-nf nextflow pipeline (nanopore/medaka workflow)
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 * Note: This runs as a shell block
 */
process ncov2019_artic_nf_pipeline_medaka {
  tag "${task.index} - ${fastq_file}"
  input:
    file fastq_file
    val ncov_prefix
    val scheme_repo_url
    val scheme_dir
    val scheme
    val scheme_version

  output:
    path "output_fasta/*.consensus.fasta", emit: ch_fasta_ncov_results
    path "ncov_output/*.qc.csv", emit: ch_qc_csv_ncov_result
    path "output_plots/*.png", emit: ch_sample_depth_ncov_results

  shell:
  '''
  ncov_out_directory="ncov_output"
  ncov_minion_medaka_out_dir="${ncov_out_directory}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"
  ncov_qc_plots_dir="${ncov_out_directory}/qc_plots"
  output_fasta="output_fasta"
  output_plots="output_plots"

  # convert nextflow variables to Bash so that the same format is used
  fastq_file=!{fastq_file}
  ncov_prefix=!{ncov_prefix}
  scheme_repo_url=!{scheme_repo_url}
  scheme_dir=!{scheme_dir}
  scheme=!{scheme}
  scheme_version=!{scheme_version}

  # move fastq file to a specific directory so that the output files will have the filename pattern:
  # <analysis_run>_<sample_id>
  # where analysis_run is ncov_prefix, and
  # sample_id is the sample UUID (=file name of the fastq file)

  # there is only one input fastq file here as ncov is executed per sample
  sample_name=`basename ${fastq_file%.*}`
  mkdir -p ${sample_name}
  mv -f ${fastq_file} ${sample_name}

  # note: we inject our configuration into ncov to override parameters
  # note: --basecalled_fastq is the directory containing the barcodes or the fastq files
  nextflow run ${PSGA_ROOT_PATH}/ncov2019-artic-nf \
      --medaka \
      --prefix ${ncov_prefix} \
      --basecalled_fastq ${sample_name} \
      --outdir ${ncov_out_directory} \
      --schemeRepoURL ${scheme_repo_url} \
      --schemeDir ${scheme_dir} \
      --scheme ${scheme} \
      --schemeVersion ${scheme_version} \
      -c ${PSGA_ROOT_PATH}/psga/ncov-custom.config \
      -c ${PSGA_ROOT_PATH}/psga/ncov-nanopore-k8s.config

  mkdir -p ${output_fasta}
  mkdir -p ${output_plots}
  mv ${ncov_minion_medaka_out_dir}/*.fasta ${output_fasta}
  mv ${ncov_minion_medaka_out_dir}/*.png ${output_plots}
  mv ${ncov_qc_plots_dir}/*.png ${output_plots}

  # this is a code correction to the nanopore medaka workflow in order to restore the correct sample names
  # in file names and file content
  for file_to_update in `ls ${output_fasta}/${ncov_prefix}_${sample_name}*.fasta`; do
      sed -i "s/${ncov_prefix}_${sample_name}/${sample_name}/" ${file_to_update}
  done

  sed -i "s/${ncov_prefix}_${sample_name}/${sample_name}/g" ${ncov_out_directory}/${ncov_prefix}.qc.csv
  mv ${ncov_out_directory}/${ncov_prefix}.qc.csv ${ncov_out_directory}/${sample_name}.qc.csv

  for file_to_rename in `find ${output_fasta} ${output_plots} -name ${ncov_prefix}_${sample_name}*`; do
      file_dir=`dirname ${file_to_rename}`
      file_name=`basename ${file_to_rename}`
      cd ${file_dir}
      rename "s/${ncov_prefix}_${sample_name}/${sample_name}/" ${file_name}
      cd - > /dev/null
  done
  '''
}


/*
 * Store ncov2019_artic output
 */
process store_ncov2019_artic_nf_output {
  publishDir "${PSGA_OUTPUT_PATH}/ncov2019-artic", mode: 'copy', overwrite: true

  input:
    path ch_ncov_fasta
    path qc_plot_image
    path qc_csv

  output:
    path ch_ncov_fasta
    path qc_plot_image
    path qc_csv

  script:

  """
  """
}

/*
 * Merge ncov QC results into one single file
 */
process merge_ncov_qc_files {
  input:
    file input_dir

  output:
    path "${output_path}", emit: ch_ncov_qc_all_samples

  script:
    output_path = "ncov_qc.csv"

  """
  # extract the header from a lineage report
  sample_qc="\$(ls *.qc.csv | head -n 1)"
  # copy over the header only (first line)
  awk 'FNR == 1' "\${sample_qc}" > ${output_path}
  # copy over the record only from all files (second line)
  awk 'FNR == 2' *.qc.csv >> ${output_path}
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
  tag "${task.index} - ${ncov_fasta}"
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
  tag "${task.index} - ${reheadered_fasta_file}"
  publishDir "${PSGA_OUTPUT_PATH}/reheadered-fasta", mode: 'copy', overwrite: true

  input:
    path reheadered_fasta_file
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
  tag "${task.index} - ${reheadered_fasta_file}"
  publishDir "${PSGA_OUTPUT_PATH}/reheadered-fasta-qc-failed", mode: 'copy', overwrite: true

  input:
    path reheadered_fasta_file
    val sample_name

  output:
    path matching_file, emit: ch_qc_failed_fasta

  script:
    matching_file = "${sample_name}.fasta"

  """
  """
}
