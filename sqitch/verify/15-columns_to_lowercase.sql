-- Verify sample_qc

DO $$
BEGIN

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample_qc' AND column_name='pct_n_bases'
        )
    );

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample_qc' AND column_name='longest_no_n_run'
        )
    );

    ASSERT (
        SELECT NOT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample' AND column_name='pct_N_bases'
        )
    );

    ASSERT (
        SELECT NOT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample' AND column_name='longest_no_N_run'
        )
    );

END $$;;
