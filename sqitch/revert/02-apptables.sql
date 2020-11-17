-- Revert ps-bahrain-covid:02-apptables from pg

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  DROP TABLE IF EXISTS "sample_comorbidity";

  DROP TABLE IF EXISTS "sample_qc";

  DROP TABLE IF EXISTS "sample";

  DROP TABLE IF EXISTS "area";

  DROP TABLE IF EXISTS "block";

  DROP TABLE IF EXISTS "comorbidity";

  DROP TYPE "gender";

  DROP TYPE "governerate";

  DROP TYPE "hospital_admittance";

COMMIT;
