-- Verify covid-pipeline:16-update_pangolin_fields om pg
DO $$
BEGIN

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='pangolin_conflict'
        )
    );

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='pangolin_ambiguity_score'
        )
    );

END $$;
