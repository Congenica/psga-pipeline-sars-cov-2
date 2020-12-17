-- Change sample.mrn type from integer to varchar
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ALTER COLUMN "mrn" TYPE VARCHAR;

COMMIT;

