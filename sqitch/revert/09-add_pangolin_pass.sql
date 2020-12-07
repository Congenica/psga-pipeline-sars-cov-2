-- Revert ps-bahrain-covid:09-add_pangolin_pass from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "pangolin_status";

    DROP TYPE "pangolin_status";

COMMIT;
