-- Verify ps-bahrain-covid:09-add_pangolin_pass on pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');

ROLLBACK;
