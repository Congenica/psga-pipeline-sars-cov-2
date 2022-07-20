/*
 * Generate some dummy output
 */
process generate_dummy_output {
  publishDir "${PSGA_OUTPUT_PATH}/merged_output", mode: 'copy', overwrite: true, pattern: 'pipeline_output.csv'

  input:
    path ch_input_files

  output:
    path "pipeline_output.csv", emit: ch_merged_csv_file

  shell:
  '''
  merged_csv_file="pipeline_output.csv"
  echo "sample_id,res1,res2" > ${merged_csv_file}

  for fq in `ls *_1.fastq.gz`; do
      sample_id=$(echo "${fq}" | cut -d'_' -f1)
      # returns the seconds and current nanoseconds since the epoch.
      res2=$(date +%s%N)
      echo "${sample_id},fastq,${res2}" >> ${merged_csv_file}
  done
  '''
}

