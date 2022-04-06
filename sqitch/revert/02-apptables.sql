-- Revert psga:02-apptables from pg

BEGIN;

  SET LOCAL search_path = psga;

  DROP TABLE IF EXISTS "analysis_run" CASCADE;

  DROP TABLE IF EXISTS "sample_qc" CASCADE;

  DROP TABLE IF EXISTS "sample" CASCADE;

COMMIT;
