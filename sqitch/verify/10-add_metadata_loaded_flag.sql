-- Verify ps-bahrain-covid:10-add_metadata_loaded_flag on pg

DO $$
BEGIN

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='metadata_loaded'
        )
    );

END $$;
