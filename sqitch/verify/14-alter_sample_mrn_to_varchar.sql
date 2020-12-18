-- Verify sample.mrn

DO $$
BEGIN
    SET LOCAL search_path = sars_cov_2;

    ASSERT (SELECT has_table_privilege('sample', 'SELECT, INSERT, UPDATE, DELETE, TRUNCATE, REFERENCES'));
    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample' AND column_name='mrn'
        )
    );

END $$;
