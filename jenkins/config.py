from jenkins.loading import load_data_from_csv


data_config = {
    # tool which generated the output file to validate
    "ncov2019_artic_nf": {
        "load_data_from_csv": load_data_from_csv,
        "config": {
            # name of the column used for listing the samples
            "sample_name_column": "sample_name",
            # map: columns in CSV vs columns in DB
            "columns_to_validate": {
                "sample_name": "sample_name",
                "pct_N_bases": "pct_n_bases",
                "pct_covered_bases": "pct_covered_bases",
                "longest_no_N_run": "longest_no_n_run",
                "num_aligned_reads": "num_aligned_reads",
                "qc_pass": "qc_pass",
            },
            # abs tolerances for the columns to round
            "columns_to_round": {
                "pct_n_bases": 0.5,
                "pct_covered_bases": 0.5,
                "longest_no_n_run": 25,
                "num_aligned_reads": 200,
            },
        },
    },
    "pangolin": {
        "load_data_from_csv": load_data_from_csv,
        "config": {
            "sample_name_column": "sample_name",
            "columns_to_validate": {
                "taxon": "sample_name",
                "lineage": "pangolin_lineage",
                "conflict": "pangolin_conflict",
                "ambiguity_score": "pangolin_ambiguity_score",
                "scorpio_call": "scorpio_call",
                "scorpio_support": "scorpio_support",
                "scorpio_conflict": "scorpio_conflict",
                "scorpio_notes": "scorpio_notes",
                "is_designated": "is_designated",
                "qc_status": "pangolin_status",
                "qc_notes": "qc_notes",
                "note": "note",
            },
            # abs tolerances for the columns to round
            "columns_to_round": {"pangolin_conflict": 0.1, "scorpio_support": 0.05, "scorpio_conflict": 0.05},
        },
    },
}
