-- Verify covid-pipeline:08-add_area_colours on pg

DO $$
DECLARE
  res int;
BEGIN
  SET LOCAL search_path = sars_cov_2;
  res := (SELECT COUNT(*) FROM area WHERE colour IS NOT NULL);
  ASSERT res = 176;
END $$;
