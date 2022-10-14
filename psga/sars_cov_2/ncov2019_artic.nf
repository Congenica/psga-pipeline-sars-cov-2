/*
 * Run: ncov2019-artic-nf nextflow pipeline (illumina workflow)
 * see: https://github.com/Congenica/ncov2019-artic-nf
 */
process ncov2019_artic_nf_pipeline_illumina {
  publishDir "${params.output_path}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: 'output_{bam,fasta,plots,variants}/*'

  tag "${task.index} - ${fastq_file}"
  input:
    path fastq_file
    val scheme_dir
    val scheme
    val scheme_version

  output:
    // retain the qc csv intentionally
    tuple path("ncov_output/*.qc.csv"), path("output_fasta/*.consensus.fa"), emit: ch_ncov_sample_fasta
    path "ncov_output/*.qc.csv", emit: ch_ncov_qc_csv
    path "output_bam/*.bam*"
    path "output_fasta/*.fa"
    path "output_plots/*.png"
    path "output_variants/*.tsv"

  shell:
  '''
  # set nextflow env vars
  export NXF_ANSI_LOG="false"

  ncov_out_dir="ncov_output"
  ncov_bam_out_dir="${ncov_out_dir}/ncovIllumina_sequenceAnalysis_trimPrimerSequences"
  ncov_fasta_out_dir="${ncov_out_dir}/ncovIllumina_sequenceAnalysis_makeConsensus"
  ncov_qc_plots_dir="${ncov_out_dir}/qc_plots"
  ncov_tsv_out_dir="${ncov_out_dir}/ncovIllumina_sequenceAnalysis_callVariants"
  bam_file_ext="mapped.primertrimmed.sorted.bam"
  output_bam="output_bam"
  output_fasta="output_fasta"
  output_plots="output_plots"
  output_variants="output_variants"

  # convert nextflow variables to Bash for convenience
  ncov_prefix=!{params.run}

  # note: we inject our configuration into ncov to override parameters
  # note: `pwd` is the workdir for this nextflow process
  # note: use /tmp as work dir, so that the intermediate files for this pipeline remain local
  #       instead of being shared with our pipeline
  nextflow run /ncov2019-artic-nf \
      --illumina \
      --prefix ${ncov_prefix} \
      --directory `eval pwd` \
      --outdir ${ncov_out_dir} \
      --schemeDir !{scheme_dir} \
      --scheme !{scheme} \
      --schemeVersion !{scheme_version} \
      -work-dir /tmp \
      -c /ncov-illumina.config

  mkdir -p ${output_bam} ${output_fasta} ${output_plots} ${output_variants}
  mv ${ncov_bam_out_dir}/*.${bam_file_ext} ${output_bam}
  mv ${ncov_fasta_out_dir}/*.primertrimmed.consensus.fa ${output_fasta}
  mv ${ncov_qc_plots_dir}/*.png ${output_plots}
  mv ${ncov_tsv_out_dir}/*.variants.tsv ${output_variants}

  # extract the sample name from one of the reads and rename the qc csv
  sample_id="$(ls *.gz | head -n 1 | cut -d"_" -f1)"
  mv ${ncov_out_dir}/${ncov_prefix}.qc.csv ${ncov_out_dir}/${sample_id}.qc.csv

  # index bam file
  samtools index ${output_bam}/${sample_id}.${bam_file_ext} ${output_bam}/${sample_id}.${bam_file_ext}.bai
  '''
}


/*
 * Run: ncov2019-artic-nf nextflow pipeline (nanopore/medaka workflow)
 * see: https://github.com/Congenica/ncov2019-artic-nf
 * Note: This runs as a shell block
 */
process ncov2019_artic_nf_pipeline_medaka {
  publishDir "${params.output_path}/ncov2019-artic", mode: 'copy', overwrite: true, pattern: 'output_{bam,fasta,plots,variants}/*'

  tag "${task.index} - ${fastq_file}"
  input:
    path fastq_file
    val scheme_dir
    val scheme
    val scheme_version

  output:
    // retain the qc csv intentionally
    tuple path("ncov_output/*.qc.csv"), path("output_fasta/*.consensus.fa"), emit: ch_ncov_sample_fasta
    path "ncov_output/*.qc.csv", emit: ch_ncov_qc_csv
    path "output_bam/*.bam*"
    path "output_fasta/*.fa"
    path "output_plots/*.png"
    path "output_variants/*.vcf.gz*"

  shell:
  '''
  # set nextflow env vars
  export NXF_ANSI_LOG="false"

  ncov_out_dir="ncov_output"
  ncov_minion_out_dir="${ncov_out_dir}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"
  ncov_qc_plots_dir="${ncov_out_dir}/qc_plots"
  bam_file_ext="primertrimmed.rg.sorted.bam"
  vcf_file_ext="pass.vcf.gz"
  output_bam="output_bam"
  output_fasta="output_fasta"
  output_plots="output_plots"
  output_variants="output_variants"

  # convert nextflow variables to Bash for convenience
  fastq_file=!{fastq_file}
  ncov_prefix=!{params.run}

  # ncov/ONT expects decompressed fastq in artic 1.1.3
  if [ "${fastq_file##*.}" = "gz" ]; then
      gunzip -c ${fastq_file} > ${fastq_file%.*}
      fastq_file="${fastq_file%.*}"
  fi

  # move fastq file to a specific directory so that the output files will have the filename pattern:
  # <analysis_run>_<sample_id>
  # where analysis_run is ncov_prefix, and
  # sample_id is the sample UUID (=file name of the fastq file)

  # there is only one input fastq file here as ncov is executed per sample
  sample_id=$( echo ${fastq_file} | cut -d '.' -f1)
  mkdir -p ${sample_id}
  mv -f ${fastq_file} ${sample_id}

  # note: we inject our configuration into ncov to override parameters
  # note: --basecalled_fastq is the directory containing the barcodes or the fastq files
  # note: use /tmp as work dir, so that the intermediate files for this pipeline remain local
  #       instead of being shared with our pipeline
  nextflow run /ncov2019-artic-nf \
      --medaka \
      --medakaModel !{params.medaka_model} \
      --normalise !{params.normalise} \
      --min_len !{params.min_len} \
      --max_len !{params.max_len} \
      --prefix ${ncov_prefix} \
      --basecalled_fastq ${sample_id} \
      --outdir ${ncov_out_dir} \
      --schemeDir !{scheme_dir} \
      --scheme !{scheme} \
      --schemeVersion !{scheme_version} \
      -work-dir /tmp \
      -c /ncov-nanopore.config

  mkdir -p ${output_bam} ${output_fasta} ${output_plots} ${output_variants}
  mv ${ncov_minion_out_dir}/*.${bam_file_ext} ${output_bam}
  mv ${ncov_minion_out_dir}/*.fasta ${output_fasta}
  mv ${ncov_qc_plots_dir}/*.png ${output_plots}
  mv ${ncov_minion_out_dir}/*.${vcf_file_ext} ${output_variants}

  # this is a code correction to the nanopore medaka workflow in order to restore the correct sample names
  # in file content and file name
  # UPDATE sample name in file content
  for file_to_update in `ls ${output_bam}/${ncov_prefix}_${sample_id}*.bam`; do
      samtools view -H ${file_to_update}  | sed "s/${ncov_prefix}_${sample_id}/${sample_id}/g" | samtools reheader - ${file_to_update} > ${file_to_update}.new
      mv ${file_to_update}.new ${file_to_update}
      samtools index ${file_to_update} ${file_to_update}.bai
  done

  for file_to_update in `ls ${output_fasta}/${ncov_prefix}_${sample_id}*.fasta`; do
      sed -i "s/${ncov_prefix}_${sample_id}/${sample_id}/" ${file_to_update}
      # use .fa extension so that this is the same as for the ncov illumina workflow
      mv ${file_to_update} ${file_to_update%.*}.fa
  done

  for file_to_update in `ls ${output_variants}/${ncov_prefix}_${sample_id}*.vcf*`; do
      zcat ${file_to_update} | sed "s/${ncov_prefix}_${sample_id}/${sample_id}/" | bgzip > "${file_to_update}.new"
      mv ${file_to_update}.new ${file_to_update}
      bcftools index -t "${file_to_update}"
  done

  sed -i "s/${ncov_prefix}_${sample_id}/${sample_id}/g" ${ncov_out_dir}/${ncov_prefix}.qc.csv
  mv ${ncov_out_dir}/${ncov_prefix}.qc.csv ${ncov_out_dir}/${sample_id}.qc.csv

  # UPDATE sample name in file name
  for file_to_rename in `find ${output_bam} ${output_fasta} ${output_plots} ${output_variants} -name ${ncov_prefix}_${sample_id}*`; do
      file_dir=`dirname ${file_to_rename}`
      file_name=`basename ${file_to_rename}`
      cd ${file_dir}
      rename "s/${ncov_prefix}_${sample_id}/${sample_id}/" ${file_name}
      cd - > /dev/null
  done
  '''
}
