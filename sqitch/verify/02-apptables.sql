-- Verify psga:02-apptables on pg

DO $$
declare
  column_count integer;
BEGIN

  SET LOCAL search_path = psga;

    -- analysis_run table verifications
  ASSERT (SELECT has_table_privilege('analysis_run', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='psga' AND table_name='analysis_run' AND column_name IN (
    'analysis_run_id',
    'analysis_run_name',
    'primer_scheme_name',
    'primer_scheme_version',
    'input_file_type',
    'workflow',
    'pipeline_version',
    'pangolin_version',
    'pangolin_data_version',
    'constellation_version',
    'scorpio_version'
  );
  ASSERT column_count = 11;

    -- sample table verifications
  ASSERT (SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='psga' AND table_name='sample' AND column_name IN (
    'sample_id',
    'analysis_run_id',
    'sample_name',
    'date_collected',
    'pangolin_lineage',
    'gisaid_id',
    'genome_length',
    'pangolin_status',
    'metadata_loaded',
    'genbank_submit_id',
    'pangolin_conflict',
    'pangolin_ambiguity_score',
    'scorpio_call',
    'scorpio_support',
    'scorpio_conflict',
    'scorpio_notes',
    'is_designated',
    'qc_notes',
    'note'
  );
  ASSERT column_count = 19;

    -- sample_qc table verifications
  ASSERT (SELECT has_table_privilege('sample_qc', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='psga' AND table_name='sample_qc' AND column_name IN (
    'sample_id',
    'pct_n_bases',
    'pct_covered_bases',
    'longest_no_n_run',
    'num_aligned_reads',
    'qc_pass',
    'qc_plot'
  );
  ASSERT column_count = 7;

END $$;
