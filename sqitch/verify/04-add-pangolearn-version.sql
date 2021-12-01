-- Verify covid-pipeline:04-add-pangolearn-version on pg

DO $$
BEGIN

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample_qc'
            AND column_name='pangolearn_version'
        )
    );

END $$;
