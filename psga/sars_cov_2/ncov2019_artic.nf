/*
 * Run: ncov2019-artic-nf nextflow pipeline (illumina workflow)
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 */
process ncov2019_artic_nf_pipeline_illumina {
  publishDir "${PSGA_OUTPUT_PATH}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: 'output_{fasta,plots}/*'

  tag "${task.index} - ${fastq_file}"
  input:
    path fastq_file
    val ncov_prefix
    val scheme_repo_url
    val scheme_dir
    val scheme
    val kit

  output:
    // retain the qc csv intentionally
    tuple path("ncov_output/*.qc.csv"), path("output_fasta/*.consensus.fa"), emit: ch_ncov_sample_fasta
    path "ncov_output/*.qc.csv", emit: ch_ncov_qc_csv
    path "output_fasta/*.fa"
    path "output_plots/*.png"

  shell:
  '''
  ncov_out_dir="ncov_output"
  ncov_fasta_out_dir="${ncov_out_dir}/ncovIllumina_sequenceAnalysis_makeConsensus"
  ncov_qc_plots_dir="${ncov_out_dir}/qc_plots"
  output_fasta="output_fasta"
  output_plots="output_plots"

  # convert nextflow variables to Bash so that the same format is used
  ncov_prefix=!{ncov_prefix}
  scheme_repo_url=!{scheme_repo_url}
  scheme_dir=!{scheme_dir}
  scheme=!{scheme}
  kit=!{kit}

  # note: we inject our configuration into ncov to override parameters
  # note: `pwd` is the workdir for this nextflow process
  # note: use /tmp as work dir, so that the intermediate files for this pipeline remain local
  #       instead of being shared with our pipeline
  nextflow run /ncov2019-artic-nf \
      --illumina \
      --prefix ${ncov_prefix} \
      --directory `eval pwd` \
      --outdir ${ncov_out_dir} \
      --schemeRepoURL ${scheme_repo_url} \
      --schemeDir ${scheme_dir} \
      --scheme ${scheme} \
      --schemeVersion ${kit} \
      -work-dir /tmp \
      -c /ncov-illumina.config

  mkdir -p ${output_fasta}
  mkdir -p ${output_plots}
  mv ${ncov_fasta_out_dir}/*.fa ${output_fasta}
  mv ${ncov_qc_plots_dir}/*.png ${output_plots}

  # extract the sample name from one of the reads and rename the qc csv
  sample_name="$(ls *.gz | head -n 1 | cut -d"_" -f1)"
  mv ${ncov_out_dir}/${ncov_prefix}.qc.csv ${ncov_out_dir}/${sample_name}.qc.csv
  '''
}


/*
 * Run: ncov2019-artic-nf nextflow pipeline (nanopore/medaka workflow)
 * see: https://github.com/connor-lab/ncov2019-artic-nf
 * Note: This runs as a shell block
 */
process ncov2019_artic_nf_pipeline_medaka {
  publishDir "${PSGA_OUTPUT_PATH}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: 'output_{fasta,plots}/*'

  tag "${task.index} - ${fastq_file}"
  input:
    path fastq_file
    val ncov_prefix
    val scheme_repo_url
    val scheme_dir
    val scheme
    val kit

  output:
    // retain the qc csv intentionally
    tuple path("ncov_output/*.qc.csv"), path("output_fasta/*.consensus.fa"), emit: ch_ncov_sample_fasta
    path "ncov_output/*.qc.csv", emit: ch_ncov_qc_csv
    path "output_fasta/*.fa"
    path "output_plots/*.png"

  shell:
  '''
  ncov_out_dir="ncov_output"
  ncov_fasta_out_dir="${ncov_out_dir}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"
  ncov_qc_plots_dir="${ncov_out_dir}/qc_plots"
  output_fasta="output_fasta"
  output_plots="output_plots"

  # convert nextflow variables to Bash so that the same format is used
  fastq_file=!{fastq_file}
  ncov_prefix=!{ncov_prefix}
  scheme_repo_url=!{scheme_repo_url}
  scheme_dir=!{scheme_dir}
  scheme=!{scheme}
  kit=!{kit}

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
  # note: use /tmp as work dir, so that the intermediate files for this pipeline remain local
  #       instead of being shared with our pipeline
  nextflow run /ncov2019-artic-nf \
      --medaka \
      --prefix ${ncov_prefix} \
      --basecalled_fastq ${sample_name} \
      --outdir ${ncov_out_dir} \
      --schemeRepoURL ${scheme_repo_url} \
      --schemeDir ${scheme_dir} \
      --scheme ${scheme} \
      --schemeVersion ${kit} \
      -work-dir /tmp \
      -c /ncov-nanopore.config

  mkdir -p ${output_fasta}
  mkdir -p ${output_plots}
  mv ${ncov_fasta_out_dir}/*.fasta ${output_fasta}
  mv ${ncov_fasta_out_dir}/*.png ${output_plots}
  mv ${ncov_qc_plots_dir}/*.png ${output_plots}

  # this is a code correction to the nanopore medaka workflow in order to restore the correct sample names
  # in file names and file content
  for file_to_update in `ls ${output_fasta}/${ncov_prefix}_${sample_name}*.fasta`; do
      sed -i "s/${ncov_prefix}_${sample_name}/${sample_name}/" ${file_to_update}
      # use .fa extension so that this is the same as for the ncov illumina workflow
      mv ${file_to_update} ${file_to_update%.*}.fa
  done

  sed -i "s/${ncov_prefix}_${sample_name}/${sample_name}/g" ${ncov_out_dir}/${ncov_prefix}.qc.csv
  mv ${ncov_out_dir}/${ncov_prefix}.qc.csv ${ncov_out_dir}/${sample_name}.qc.csv

  for file_to_rename in `find ${output_fasta} ${output_plots} -name ${ncov_prefix}_${sample_name}*`; do
      file_dir=`dirname ${file_to_rename}`
      file_name=`basename ${file_to_rename}`
      cd ${file_dir}
      rename "s/${ncov_prefix}_${sample_name}/${sample_name}/" ${file_name}
      cd - > /dev/null
  done
  '''
}
