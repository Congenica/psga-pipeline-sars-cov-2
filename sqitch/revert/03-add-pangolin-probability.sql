-- Revert ps-bahrain-covid:03-add-pangolin-probability from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "pangolin_probability";

COMMIT;
