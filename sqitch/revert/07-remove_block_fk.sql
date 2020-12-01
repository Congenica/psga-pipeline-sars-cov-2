-- Revert ps-bahrain-covid:07-remove_block_fk from pg

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  ALTER TABLE sample ADD CONSTRAINT
      sample_block_number_fkey FOREIGN KEY (block_number) REFERENCES block(number) DEFERRABLE INITIALLY DEFERRED;

  COMMENT ON COLUMN "sample"."block_number" IS
    'foreign key to block table, low level geolocation from where sample was collected';

COMMIT;
