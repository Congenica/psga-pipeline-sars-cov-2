from sqlalchemy.orm import scoped_session

from scripts.db.models import AnalysisRun, Sample


def get_analysis_run(session: scoped_session, analysis_run_name: str) -> AnalysisRun:
    """
    Return the analysis run with this name or None
    """
    analysis_run = (
        session.query(AnalysisRun)
        .filter(
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )
    return analysis_run


def get_analysis_run_sample(session: scoped_session, analysis_run_name: str, sample_name: str) -> AnalysisRun:
    """
    Return the sample with this name for this analysis run name or None
    """
    sample = (
        session.query(Sample)
        .join(AnalysisRun)
        .filter(
            Sample.sample_name == sample_name,
            AnalysisRun.analysis_run_name == analysis_run_name,
        )
        .one_or_none()
    )
    return sample
