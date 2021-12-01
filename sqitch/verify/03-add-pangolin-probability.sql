-- Verify covid-pipeline:03-add-pangolin-probability om pg
DO $$
BEGIN

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='pangolin_probability'
        )
    );

END $$;
