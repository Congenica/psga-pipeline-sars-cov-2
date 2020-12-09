-- Deploy ps-bahrain-covid:09-add_pangolin_pass to pg
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    CREATE TYPE "pangolin_status" AS ENUM (
        'unknown',
        'fail'
       ,'passed_qc'
      );

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "pangolin_status" "pangolin_status" NOT NULL DEFAULT 'unknown';

    COMMENT ON COLUMN "sample"."pangolin_status"
        IS 'Reported pangolin lineage status';

COMMIT;
