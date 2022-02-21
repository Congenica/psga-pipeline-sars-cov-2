-- Revert covid-pipeline:02-apptables from pg

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  DROP TABLE IF EXISTS "analysis_run";

  DROP TABLE IF EXISTS "sample_qc";

  DROP TABLE IF EXISTS "sample";

COMMIT;
