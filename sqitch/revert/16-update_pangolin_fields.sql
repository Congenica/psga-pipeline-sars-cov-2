-- Revert covid-pipeline:16-update_pangolin_fields from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "pangolin_conflict";
    ALTER TABLE "sample" DROP COLUMN IF EXISTS "pangolin_ambiguity_score";

COMMIT;

