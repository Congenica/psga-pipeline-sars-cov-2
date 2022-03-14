-- Verify psga:appschema on pg

BEGIN;

    SELECT pg_catalog.has_schema_privilege('psga', 'usage');

ROLLBACK;
