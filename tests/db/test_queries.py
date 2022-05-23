from scripts.db.queries import get_analysis_run, get_analysis_run_sample


def test_get_analysis_run(db_session, populated_db_session_with_sample):
    expected_analysis_run_name = "just_a_name"

    analysis_run = get_analysis_run(db_session, "fake_analysis_run")
    assert analysis_run is None

    analysis_run = get_analysis_run(db_session, "just_a_name")

    assert analysis_run is not None
    assert analysis_run.analysis_run_name == expected_analysis_run_name


def test_get_analysis_run_sample(db_session, populated_db_session_with_sample):
    expected_analysis_run_name = "just_a_name"
    expected_analysis_run_sample_name = "7284954"

    sample = get_analysis_run_sample(db_session, "fake_analysis_run", "fake_sample")
    assert sample is None

    sample = get_analysis_run_sample(db_session, "fake_analysis_run", expected_analysis_run_sample_name)
    assert sample is None

    sample = get_analysis_run_sample(db_session, expected_analysis_run_name, "fake_sample")
    assert sample is None

    sample = get_analysis_run_sample(db_session, expected_analysis_run_name, expected_analysis_run_sample_name)
    assert sample is not None
    assert sample.sample_name == expected_analysis_run_sample_name
