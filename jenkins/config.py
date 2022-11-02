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
                "PRIMER_QC",
                "PRIMER_NUMREADS",
                "PRIMER_COVBASES",
                "PRIMER_COVERAGE",
                "PCT_N_BASES",
                "PCT_COVERED_BASES",
                "LONGEST_NO_N_RUN",
                "QC_PASS",
                "LINEAGE",
                "CONFLICT",
                "AMBIGUITY_SCORE",
                "SCORPIO_CALL",
                "SCORPIO_SUPPORT",
                "SCORPIO_CONFLICT",
                "SCORPIO_NOTES",
                "PANGOLIN_VERSION",
                "SCORPIO_VERSION",
                "CONSTELLATION_VERSION",
                "IS_DESIGNATED",
                "QC_STATUS",
                "QC_NOTES",
            ],
            # abs tolerances for the columns to round
            "columns_to_round": {
                "PCT_N_BASES": 0.5,
                "PCT_COVERED_BASES": 0.5,
                "LONGEST_NO_N_RUN": 25,
                "CONFLICT": 0.1,
                "SCORPIO_SUPPORT": 0.05,
                "SCORPIO_CONFLICT": 0.05,
            },
        },
    },
}
