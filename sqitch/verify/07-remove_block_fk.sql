-- Verify covid-pipeline:07-remove_block_fk on pg

DO $$
DECLARE
  res int;
BEGIN
  res := (SELECT COUNT(*) FROM pg_catalog.pg_constraint WHERE conname = 'sample_block_number_fkey');
  ASSERT res = 0;
END $$;
