data_config = {
    # pathogen
    "sars_cov_2": {
        "config": {
            # name of the column used for listing the samples
            "sample_name_column": "SAMPLE_ID",
            "columns_to_validate": [
                "SAMPLE_ID",
                "STATUS",
                "CONTAMINATED_READS",
                "PRIMER_INPUT",
                "PRIMER_DETECTED",
                "PRIMER_NUMREADS",
                "PRIMER_UNIQUE_NUMREADS",
                "PRIMER_AMBIGUOUS_NUMREADS",
                "PRIMER_COVERAGE",
                "NCOV_PCT_N_BASES",
                "NCOV_PCT_COVERED_BASES",
                "NCOV_LONGEST_NO_N_RUN",
                "NCOV_QC_PASS",
                "PANGOLIN_LINEAGE",
                "PANGOLIN_CONFLICT",
                "PANGOLIN_AMBIGUITY_SCORE",
                "PANGOLIN_SCORPIO_CALL",
                "PANGOLIN_SCORPIO_SUPPORT",
                "PANGOLIN_SCORPIO_CONFLICT",
                "PANGOLIN_SCORPIO_NOTES",
                "PANGOLIN_VERSION",
                "PANGOLIN_SCORPIO_VERSION",
                "PANGOLIN_CONSTELLATION_VERSION",
                "PANGOLIN_IS_DESIGNATED",
                "PANGOLIN_QC_STATUS",
                "PANGOLIN_QC_NOTES",
            ],
            # abs tolerances for the columns to round
            "columns_to_round": {
                "NCOV_PCT_N_BASES": 0.5,
                "NCOV_PCT_COVERED_BASES": 0.5,
                "NCOV_LONGEST_NO_N_RUN": 25,
                "PANGOLIN_CONFLICT": 0.1,
                "PANGOLIN_SCORPIO_SUPPORT": 0.05,
                "PANGOLIN_SCORPIO_CONFLICT": 0.05,
            },
        },
    },
}
