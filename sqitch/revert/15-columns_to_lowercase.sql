-- revert column names to original names
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" RENAME "pct_n_bases" TO "pct_N_bases";

    ALTER TABLE "sample_qc" RENAME "longest_no_n_run" TO "longest_no_N_run";

COMMIT;
