-- Deploy covid-pipeline:10-add_metadata_loaded_flag to pg
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "metadata_loaded" BOOLEAN;

    COMMENT ON COLUMN "sample"."metadata_loaded"
        IS 'Metadata was loaded from Bahrain medical records system I-SEHA';

COMMIT;
