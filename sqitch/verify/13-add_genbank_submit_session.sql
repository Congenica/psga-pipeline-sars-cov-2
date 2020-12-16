-- Verify ps-bahrain-covid:13-add_genbank_submit_session om pg

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    SELECT genbank_submit_id FROM sample WHERE FALSE;

    SELECT 1/(COUNT(*)-1) FROM information_schema.columns WHERE table_name='sample' and column_name='genbank_id';

COMMIT;

