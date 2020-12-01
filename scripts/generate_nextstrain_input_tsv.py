from typing import List
from datetime import datetime
from dataclasses import dataclass
import csv

import click
from sqlalchemy.orm import joinedload

from db.database import session_handler
from db.models import Sample, SampleQC


@dataclass
class NextstrainSampleMetadataInput:
    strain: str
    date: datetime
    division: str
    area: str
    country_exposure: str
    division_exposure: str
    length: int
    age: int
    sex: str
    pangolin_lineage: str
    nationality: str
    ct_value: float
    symptoms: str
    comorbitities: str
    travel_exposure: str
    hospital_admittance: str

    gisaid_epi_isl = None
    genbank_accession = None
    region_exposure = None
    authors = None
    title = None
    date_submitted = None
    virus: str = "SARS-CoV-2"
    region: str = "Asia"
    country: str = "Bahrain"
    segment: str = "genome"
    host: str = "human"
    originating_lab: str = "Bahrain Public Health Laboratory"
    submitting_lab: str = "Bahrain Public Health Laboratory"


def get_data_for_nextstrain() -> List[NextstrainSampleMetadataInput]:
    res = []
    with session_handler() as session:
        # We process only the samples, which were marked by ncov pipeline as qc_pass=TRUE
        samples = (
            session.query(Sample)
            .options(joinedload(Sample.comorbidities))
            .join(Sample.sample_qc)
            .filter(SampleQC.qc_pass)
            .all()
        )
        for sample in samples:
            division = sample.governorate_name
            country_exposure = sample.travel_exposure
            if not country_exposure:
                country_exposure = NextstrainSampleMetadataInput.country
            res.append(
                NextstrainSampleMetadataInput(
                    strain=sample.lab_id,
                    date=sample.date_collected,
                    division=division,
                    area=sample.area_name,
                    country_exposure=country_exposure,
                    division_exposure=division,
                    length=sample.genome_length,
                    age=sample.age,
                    sex=sample.gender,
                    pangolin_lineage=sample.pangolin_lineage,
                    nationality=sample.nationality,
                    ct_value=sample.ct_value,
                    symptoms=sample.symptoms,
                    comorbitities=",".join([x.description for x in sample.comorbidities]),
                    travel_exposure=sample.travel_exposure,
                    hospital_admittance=sample.hospital_admittance,
                )
            )
    return res


@click.command()
@click.option(
    "--output",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="output tsv file path, where sample data will be written to",
)
def generate_nextstrain_input_tsv(output: str) -> None:
    """
    Generate a tsv file from a database, which will be used in nextstrain pipeline
    """
    samples = get_data_for_nextstrain()

    with open(output, "w") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        header_keys = None
        for sample in samples:
            if header_keys is None:
                header_keys = sample.__dict__.keys()
                tsv_writer.writerow(header_keys)
            tsv_writer.writerow([sample.__getattribute__(key) for key in header_keys])


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_nextstrain_input_tsv()
