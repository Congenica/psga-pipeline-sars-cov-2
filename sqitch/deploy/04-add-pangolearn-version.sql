-- Add pangolearn_version to sample_qc table
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" ADD COLUMN IF NOT EXISTS "pangolearn_version" VARCHAR;

    COMMENT ON COLUMN "sample_qc"."pangolearn_version"
        IS 'The version of pangoLEARN used by Pangolin.';
COMMIT;
