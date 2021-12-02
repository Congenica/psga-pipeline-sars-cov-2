-- Deploy covid-pipeline:05-governorates to pg
-- requires: 02-apptables

BEGIN;

  SET LOCAL search_path = sars_cov_2;

  INSERT INTO governorate VALUES ('Capital'::governorate_name, 26.23269, 50.57811);
  INSERT INTO governorate VALUES ('Northern'::governorate_name, 26.096741, 50.500927);
  INSERT INTO governorate VALUES ('Southern'::governorate_name, 25.99, 50.5575);
  INSERT INTO governorate VALUES ('Muharraq'::governorate_name, 26.25278, 50.63632);

COMMIT;
