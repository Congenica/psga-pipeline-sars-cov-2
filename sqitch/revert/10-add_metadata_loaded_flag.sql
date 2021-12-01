-- Revert covid-pipeline:10-add_metadata_loaded_flag from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "metadata_loaded";

COMMIT;
