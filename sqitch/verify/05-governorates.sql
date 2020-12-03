-- Verify ps-bahrain-covid:05-governorates on pg

DO $$
DECLARE
  res int;
BEGIN
  SET LOCAL search_path = sars_cov_2;
  res := (SELECT COUNT(*) FROM governorate);
  ASSERT res = 4;
END $$;
