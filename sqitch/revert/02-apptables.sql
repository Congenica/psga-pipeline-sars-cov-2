-- Revert psga:02-apptables from pg

BEGIN;

  SET LOCAL search_path = psga;

  DROP TABLE IF EXISTS "analysis_run";

  DROP TABLE IF EXISTS "sample_qc";

  DROP TABLE IF EXISTS "sample";

COMMIT;
