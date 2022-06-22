data_config = {
    # pathogen
    "sars_cov_2": {
        "config": {
            # name of the column used for listing the samples
            "sample_name_column": "sample_id",
            "columns_to_validate": [
                "sample_id",
                "pct_N_bases",
                "pct_covered_bases",
                "longest_no_N_run",
                "num_aligned_reads",
                "qc_pass",
                "lineage",
                "conflict",
                "ambiguity_score",
                "scorpio_call",
                "scorpio_support",
                "scorpio_conflict",
                "scorpio_notes",
                "is_designated",
                "qc_status",
                "qc_notes",
                "note",
            ],
            # abs tolerances for the columns to round
            "columns_to_round": {
                "pct_N_bases": 0.5,
                "pct_covered_bases": 0.5,
                "longest_no_N_run": 25,
                "num_aligned_reads": 200,
                "conflict": 0.1,
                "scorpio_support": 0.05,
                "scorpio_conflict": 0.05,
            },
        },
    },
}
