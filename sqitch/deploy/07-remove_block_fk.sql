-- Deploy ps-bahrain-covid:07-remove_block_fk to pg
-- requires: 02-apptables

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  ALTER TABLE sample DROP CONSTRAINT sample_block_number_fkey;

  COMMENT ON COLUMN "sample"."block_number" IS
    'low level geolocation from where sample was collected';

COMMIT;
