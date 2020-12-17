-- Change column names to lowercase
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample_qc" RENAME "pct_N_bases" TO "pct_n_bases";

    ALTER TABLE "sample_qc" RENAME "longest_no_N_run" TO "longest_no_n_run";

COMMIT;
