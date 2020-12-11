-- Revert ps-bahrain-covid:12-add-nucleotide-muts-column from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "nucleotide_muts";

COMMIT;

