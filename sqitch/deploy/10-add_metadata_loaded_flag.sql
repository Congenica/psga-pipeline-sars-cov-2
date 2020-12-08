-- Deploy ps-bahrain-covid:10-add_metadata_loaded_flag to pg
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" ADD COLUMN IF NOT EXISTS "metadata_loaded" BOOLEAN DEFAULT FALSE;

    COMMENT ON COLUMN "sample_qc"."metadata_loaded"
        IS 'Metadata was loaded from Bahrain medical records system I-SEHA';

COMMIT;
