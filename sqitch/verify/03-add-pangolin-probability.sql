-- Verify ps-bahrain-covid:03-add-pangolin-probability om pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES');

COMMIT;