/*
 * Prepare input tsv for microreact
 */
process prepare_microreact_tsv {
  input:
    path ncov_qc_to_db_submit_completion_flag
    path pangolin_to_db_submit_completion_flag

  output:
    path params.microreact_tsv

  """
  python /app/scripts/generate_microreact_input.py --output ${params.microreact_tsv}
  """
}