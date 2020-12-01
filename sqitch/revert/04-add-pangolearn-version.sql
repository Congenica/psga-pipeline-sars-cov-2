-- Revert ps-bahrain-covid:04-add-pangolearn-version from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" DROP COLUMN IF EXISTS "pangolearn_version";

COMMIT;
