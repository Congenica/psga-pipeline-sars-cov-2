-- Verify ps-bahrain-covid:06-blocks on pg

DO $$
DECLARE
  res int;
BEGIN
  SET LOCAL search_path = sars_cov_2;
  res := (SELECT COUNT(*) FROM area);
  ASSERT res = 176;
END $$;
