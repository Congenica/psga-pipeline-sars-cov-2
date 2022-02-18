-- Revert covid-pipeline:17-add_scorpio_fields from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "scorpio_call";
    ALTER TABLE "sample" DROP COLUMN IF EXISTS "scorpio_support";
    ALTER TABLE "sample" DROP COLUMN IF EXISTS "scorpio_conflict";
    ALTER TABLE "sample" DROP COLUMN IF EXISTS "note";

COMMIT;
