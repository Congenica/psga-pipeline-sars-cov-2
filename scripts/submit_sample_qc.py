import click
from bahrain_covid.database import session_handler
from bahrain_covid.models import Sample, SampleQc


@click.command()
@click.option("--sample_name", type=str, required=True, help="lab sample identifier")
@click.option("--pct_n_bases", type=float, required=True, help="percentage of N bases in sequenced genome")
@click.option("--pct_covered_bases", type=float, required=True, help="percentage of covered bases in sequenced genome")
@click.option("--longest_no_n_run", type=int, required=True, help="longest sequence run without N base")
@click.option("--num_aligned_reads", type=int, required=True, help="number of reads succesfully aligned and mapped")
@click.option("--qc_pass", type=bool, required=True, help="Has sample passed QC")
@click.option("--qc_plot", type=click.Path(exists=True, file_okay=True, readable=True), required=True,
              help="QC plot image")
@click.option("--pipeline_version", type=str, required=True, help="mapping pipeline version")
def submit_sample_qc(sample_name: str, pct_n_bases: float, pct_covered_bases: float, longest_no_n_run: int,
                     num_aligned_reads: int, qc_pass: bool, qc_plot: str, pipeline_version: str) -> None:
    """
    Submit samples QC to the database, generated by ncov pipeline
    """
    with session_handler() as session:
        sample = session.query(Sample).filter_by(lab_id=sample_name).join(Sample.sample_qc).one_or_none()
        if not sample:
            sample = Sample(lab_id=sample_name)
            session.add(sample)
        if not sample.sample_qc:
            sample.sample_qc = SampleQc()
        sample_qc = sample.sample_qc

        sample_qc.pct_N_bases = pct_n_bases
        sample_qc.pct_covered_bases = pct_covered_bases
        sample_qc.longest_no_N_run = longest_no_n_run
        sample_qc.num_aligned_reads = num_aligned_reads
        sample_qc.qc_pass = qc_pass

        with open(qc_plot, "rb") as f:
            sample_qc.qc_plot = bytearray(f.read())

        sample_qc.pipeline_version = pipeline_version


if __name__ == "__main__":
    submit_sample_qc()
