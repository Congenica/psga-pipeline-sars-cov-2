-- Verify covid-pipeline:appschema on pg

BEGIN;

    SELECT pg_catalog.has_schema_privilege('sars_cov_2', 'usage');

ROLLBACK;
