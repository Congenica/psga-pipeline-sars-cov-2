-- Revert psga:appschema from pg

BEGIN;

    DROP SCHEMA psga;

COMMIT;
