-- Add amino acid mutation column to sample table
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "amino_acid_muts" VARCHAR;

    COMMENT ON COLUMN "sample"."amino_acid_muts"
        IS 'The amino acid mutations of the sample as calculated by Nextstrain.';
COMMIT;
