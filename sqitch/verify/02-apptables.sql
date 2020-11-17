-- Verify ps-bahrain-covid:02-apptables on pg

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  SELECT has_table_privilege('area', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');
  SELECT has_table_privilege('block', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');
  SELECT has_table_privilege('comorbidity', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');
  SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');
  SELECT has_table_privilege('sample_qc', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');
  SELECT has_table_privilege('sample_comorbidity', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');

ROLLBACK;
