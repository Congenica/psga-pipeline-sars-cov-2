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
  WHERE table_schema='sars_cov_2' AND table_name='sample' AND column_name IN ('sample_id', 'lab_id', 'date_collected', 'data_sequenced', 'sequencing_run', 'sample_number', 'pangolin_lineage', 'gisaid_id', 'genome_length');
  ASSERT column_count = 20;

    -- sample_qc table verifications
  ASSERT (SELECT has_table_privilege('sample_qc', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='sample_qc' AND column_name IN ('sample_id', 'pct_covered_bases', 'num_aligned_reads', 'qc_pass', 'qc_plot', 'pipeline_version');
  ASSERT column_count = 6;

END $$;
