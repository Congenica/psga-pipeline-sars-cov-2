-- Verify ps-bahrain-covid:11-add-amino-acid-muts-column om pg

DO $$
BEGIN

    ASSERT (
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.columns
            WHERE table_schema='sars_cov_2' AND table_name='sample'
            AND column_name='amino_acid_muts'
        )
    );

END $$;
