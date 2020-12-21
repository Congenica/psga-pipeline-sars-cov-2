/*
 * Concatenate reheadered fasta files for Nextstrain pipeline
 */
process concatenate_fasta {
  input:
    val root_genome_fasta
    path reheadered_fasta
    path archived_fasta

  output:
    path output_file

  script:
    files_dir = "./"
    output_file = "nextstrain.fasta"

  """
  # create links in ${files_dir} to the archived files, so that these can be concatenated
  python /app/scripts/link_archived_fasta.py --destination ${files_dir}

  echo "FASTA files to concatenate:"
  ls -l

  python /app/scripts/concatenate_fasta.py --output ${output_file} --root-genome ${root_genome_fasta} ${files_dir}
  """
}

/*
 * Prepare metadata tsv file, which is used as an input for nextstrain pipeline
 */
process prepare_tsv_for_nextstrain {
  input:
    path ncov_qc_to_db_submit_completion_flag
    path pangolin_to_db_submit_completion_flag

  output:
    path nextstrain_analysis_tsv

  script:
    nextstrain_analysis_tsv = "nextstrain_metadata.tsv"

  """
  python /app/scripts/generate_nextstrain_input_tsv.py --output ${nextstrain_analysis_tsv}
  """
}

/*
 * Run: nextstrain-ncov snakemake pipeline
 * see: https://github.com/nextstrain/ncov
 */
process nextstrain_pipeline {
  input:
    path metadata
    path fasta

  output:
    path "${nextstrain_out_directory}/*", emit: ch_all_nextstrain_results
    path "${nextstrain_out_directory}/bahrain/ncov_with_accessions.json", emit: ch_nextstrain_ncov_with_accessions_json
    path "${nextstrain_out_directory}/bahrain/aa_muts.json", emit: ch_nextstrain_aa_muts_json
    path "${nextstrain_out_directory}/bahrain/nt_muts.json", emit: ch_nextstrain_nt_muts_json
    path "${nextstrain_out_directory}/bahrain/tree.nwk", emit: ch_nextstrain_tree_nwk

  script:
    nextstrain_out_directory = "nextstrain_output"
  """
  # The snakemake custom profile is configured so that the input data are in: /nextstrain/data/
  # These copies must be executed as root as they alter the container, hence --user 0:0 in nextflow.config.
  # note: ln -s does not work here
  cp ${metadata} /nextstrain/data/nextstrain_metadata.tsv
  cp ${fasta} /nextstrain/data/nextstrain.fasta

  # create the output dir. This must be done on /
  mkdir -p ${nextstrain_out_directory}
  cd /nextstrain
  # run the nextstrain Snakemake pipeline
  snakemake --profile /custom_profile

  # copy the output files to the caller and set user permissions. This must be done on /
  # note: ln -s does not work here
  cd -
  cp -R /nextstrain/results/* ${nextstrain_out_directory}
  chown -R \${UID}:\${GID} ${nextstrain_out_directory}
  """
}

/*
 * Store nextstrain output
 * We publish only a single channel. This way multiple channels won't conflict on publish
 */
process store_nextstrain_output {
  publishDir "${COVID_PIPELINE_NEXTSTRAIN_PATH}/${workflow.sessionId}", mode: 'copy', overwrite: true
  publishDir "${COVID_PIPELINE_NEXTSTRAIN_PATH}/latest", mode: 'copy', overwrite: true

  input:
    path ch_all_nextstrain_results

  output:
    path ch_all_nextstrain_results

  script:

  """
  """
}

/*
 * Load Nextstrain data to the database.
 * Data are: amino acid and nucleotide mutations.
 */
process load_nextstrain_data_to_db {
  input:
    file(ch_nextstrain_aa_muts_json_file)
    file(ch_nextstrain_nt_muts_json_file)
    file(ch_nextstrain_tree_nwk_file)

  output:
    path ch_load_nextstrain_data_done, emit: ch_nextstrain_data_load_done

  script:
    ch_load_nextstrain_data_done = "load_nextstrain_data_to_db.done"

  """
  # load amino acid mutations
  python /app/scripts/load_nextstrain_aa_muts_to_db.py \
    --aa-muts-json "${ch_nextstrain_aa_muts_json_file}" \
    --tree-nwk "${ch_nextstrain_tree_nwk_file}"

  # load nucleotide mutations
  python /app/scripts/load_nextstrain_nt_muts_to_db.py \
    --nt-muts-json "${ch_nextstrain_nt_muts_json_file}" \
    --tree-nwk "${ch_nextstrain_tree_nwk_file}"

  touch ${ch_load_nextstrain_data_done}
  """
}
