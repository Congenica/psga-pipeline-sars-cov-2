-- Add nucleotide mutation column to sample table
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "nucleotide_muts" VARCHAR;

    COMMENT ON COLUMN "sample"."nucleotide_muts"
        IS 'The nucleotide mutations of the sample as calculated by Nextstrain.';
COMMIT;
