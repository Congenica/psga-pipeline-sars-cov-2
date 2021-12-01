-- Verify covid-pipeline:02-apptables on pg

DO $$
declare
  column_count integer;
BEGIN

  SET LOCAL search_path = sars_cov_2;

  -- area table verifications
  ASSERT (SELECT has_table_privilege('area', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='area' AND column_name IN ('name','latitude','longitude','colour');
  ASSERT column_count = 4;

  -- block table verifications
  ASSERT (SELECT has_table_privilege('block', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='block' AND column_name IN ('number','latitude','longitude','colour');
  ASSERT column_count = 4;

    -- comorbidity table verifications
  ASSERT (SELECT has_table_privilege('comorbidity', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='comorbidity' AND column_name IN ('comorbidity_id','description');
  ASSERT column_count = 2;

    -- sample table verifications
  ASSERT (SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='sample' AND column_name IN ('sample_id', 'lab_id', 'date_collected', 'data_sequenced', 'sequencing_run', 'gender', 'age', 'nationality', 'governorate_name', 'area_name', 'block_number', 'sample_number', 'ct_value', 'symptoms', 'travel_exposure', 'hospital_admittance', 'pangolin_lineage', 'gisaid_id', 'genome_length', 'mrn');
  ASSERT column_count = 20;

    -- sample_qc table verifications
  ASSERT (SELECT has_table_privilege('sample_qc', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='sample_qc' AND column_name IN ('sample_id', 'pct_covered_bases', 'num_aligned_reads', 'qc_pass', 'qc_plot', 'pipeline_version');
  ASSERT column_count = 6;

    -- sample_comorbidity table verifications
  ASSERT (SELECT has_table_privilege('sample_comorbidity', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='sample_comorbidity' AND column_name IN ('sample_id','comorbidity_id');
  ASSERT column_count = 2;

    -- governorate table verifications
  ASSERT (SELECT has_table_privilege('governorate', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
  SELECT count(*) INTO column_count
  FROM information_schema.columns
  WHERE table_schema='sars_cov_2' AND table_name='governorate' AND column_name IN ('name','latitude','longitude','colour');
  ASSERT column_count = 4;

END $$;
