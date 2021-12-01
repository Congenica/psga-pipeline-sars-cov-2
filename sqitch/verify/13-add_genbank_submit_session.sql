-- Verify covid-pipeline:13-add_genbank_submit_session om pg

DO $$
BEGIN

    SET LOCAL search_path = sars_cov_2;

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='genbank_submit_id'
        )
    );

    ASSERT (
        SELECT NOT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='genbank_id'
        )
    );

END $$;
