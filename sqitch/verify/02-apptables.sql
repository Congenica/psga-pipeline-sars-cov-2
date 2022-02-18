-- Verify covid-pipeline:02-apptables on pg

DO $$
declare
  column_count integer;
BEGIN

  SET LOCAL search_path = sars_cov_2;

    -- sample table verifications
  ASSERT (SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='sample' AND column_name IN (
    'sample_id',
    'sample_name',
    'date_collected',
    'pangolin_lineage',
    'gisaid_id',
    'genome_length',
    'pangolin_status',
    'metadata_loaded',
    'genbank_submit_id',
    'pangolin_conflict',
    'pangolin_ambiguity_score'
  );
  ASSERT column_count = 11;

    -- sample_qc table verifications
  ASSERT (SELECT has_table_privilege('sample_qc', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='sample_qc' AND column_name IN (
    'sample_id',
    'pct_n_bases',
    'pct_covered_bases',
    'longest_no_n_run',
    'num_aligned_reads',
    'qc_pass',
    'qc_plot',
    'pipeline_version',
    'pangolearn_version',
    'pangolin_version',
    'pango_version'
  );
  ASSERT column_count = 11;

END $$;
