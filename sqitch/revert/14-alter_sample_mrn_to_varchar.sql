-- Revert sample.mrn type from varchar to integer

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ALTER COLUMN "mrn" TYPE INTEGER USING mrn::integer;

COMMIT;

