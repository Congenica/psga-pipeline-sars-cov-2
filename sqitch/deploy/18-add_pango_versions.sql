-- Add pango versions to sample_qc table
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" ADD COLUMN IF NOT EXISTS "pangolin_version" VARCHAR;
    ALTER TABLE "sample_qc" ADD COLUMN IF NOT EXISTS "pango_version" VARCHAR;

    COMMENT ON COLUMN "sample_qc"."pangolin_version"
        IS 'The version of Pangolin.';
    COMMENT ON COLUMN "sample_qc"."pango_version"
        IS 'The version of pango used by Pangolin.';
COMMIT;
