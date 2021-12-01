-- Revert covid-pipeline:05-governorates from pg

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  DELETE FROM governorate WHERE name in (
    'Capital'::governorate_name,
    'Northern'::governorate_name,
    'Southern'::governorate_name,
    'Muharraq'::governorate_name
  );

COMMIT;
