-- Add pangolin_probability to sample table
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "pangolin_probability" DOUBLE PRECISION;

    COMMENT ON COLUMN "sample"."pangolin_probability"
        IS 'The probability of the sample as calculated by Pangolin.';
COMMIT;
