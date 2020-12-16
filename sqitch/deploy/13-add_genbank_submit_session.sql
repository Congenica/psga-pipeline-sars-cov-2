-- Add genbank batch submission ID for the sample. Acts as an indicator, that the sample was already submitted
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "genbank_submit_id" VARCHAR;

    COMMENT ON COLUMN "sample"."genbank_submit_id"
        IS 'Unique identifier of GenBank submission session, which was used to submit a sample.';

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "genbank_id";
COMMIT;
