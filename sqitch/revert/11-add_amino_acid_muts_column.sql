-- Revert ps-bahrain-covid:11-add-amino-acid-muts-column from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "amino_acid_muts";

COMMIT;

