-- Revert covid-pipeline:04-add-pangolearn-version from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" DROP COLUMN IF EXISTS "pangolin_version";
    ALTER TABLE "sample_qc" DROP COLUMN IF EXISTS "pango_version";

COMMIT;
