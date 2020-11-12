-- Revert ps-bahrain-covid:appschema from pg

BEGIN;

    DROP SCHEMA sars_cov_2;

COMMIT;
