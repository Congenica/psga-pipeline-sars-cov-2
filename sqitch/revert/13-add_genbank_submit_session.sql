-- Revert ps-bahrain-covid:13-add_genbank_submit_session from pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "genbank_submit_id";

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "genbank_id" INTEGER;

    COMMENT ON COLUMN "sample"."genbank_id" IS 'GeneBank accession of virus genome';

COMMIT;

