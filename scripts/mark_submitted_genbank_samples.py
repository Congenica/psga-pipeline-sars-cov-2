import click

from scripts.db.database import session_handler
from scripts.db.models import Sample


@click.command()
@click.option(
    "--sample_names_txt",
    type=click.Path(file_okay=True, readable=True, exists=True),
    required=True,
    help="A text file, containing sample names in each of the lines",
)
@click.option(
    "--submit_id",
    type=str,
    required=True,
    help="Unique identifier of the submission to GenBank",
)
def mark_submitted_genbank_samples(
    sample_names_txt: str,
    submit_id: str,
) -> None:
    """
    Mark samples as submitted to genbank with submission identifier.
    Marked samples will not be submitted to GenBank in the future
    """
    with open(sample_names_txt, "r") as f:
        sample_names = f.readlines()
    sample_names = [sample.strip() for sample in sample_names]

    with session_handler() as session:
        samples = session.query(Sample).filter(Sample.sample_name.in_(sample_names)).all()

        for sample in samples:
            sample.genbank_submit_id = submit_id


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    mark_submitted_genbank_samples()
