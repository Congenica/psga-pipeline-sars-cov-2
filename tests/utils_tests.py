from scripts.db.models import AnalysisRun, Sample


def assert_files_are_equal(file1, file2):
    with open(file1, "r") as f1:
        with open(file2, "r") as f2:
            content_1 = f1.read()
            content_2 = f2.read()
            assert content_1 == content_2


def get_analysis_run_samples(db_session, analysis_run_name):
    samples = (
        db_session.query(Sample)
        .join(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .all()
    )

    return samples
