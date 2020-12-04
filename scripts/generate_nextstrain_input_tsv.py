from typing import List
from dataclasses import dataclass
import csv

import click
from sqlalchemy.orm import joinedload

from db.database import session_handler
from db.models import Sample, SampleQC, GenderEnum, HospitalAdmittanceEnum


@dataclass
class NextstrainSampleMetadataInput:
    strain: str
    date: str
    country: str
    division: str
    area: str
    country_exposure: str
    division_exposure: str
    length: int
    age: int
    sex: GenderEnum
    originating_lab: str
    submitting_lab: str
    pangolin_lineage: str
    nationality: str
    ct_value: float
    symptoms: str
    comorbitities: str
    travel_exposure: str
    hospital_admittance: HospitalAdmittanceEnum

    gisaid_epi_isl = None
    genbank_accession = None
    region_exposure = None
    authors = None
    title = None
    date_submitted = None
    virus: str = "SARS-CoV-2"
    region: str = "Asia"
    segment: str = "genome"
    host: str = "human"


root_genome = NextstrainSampleMetadataInput(
    strain="NC_045512.2",
    date="2019-12-20",
    country="China",
    division="C",
    area="",
    country_exposure="",
    division_exposure="",
    length=29903,
    age=94,
    sex=GenderEnum.F,
    originating_lab="China",
    submitting_lab="WHO",
    pangolin_lineage="A",
    nationality="China",
    ct_value=0.2657696717449,
    symptoms="",
    comorbitities="",
    travel_exposure="",
    hospital_admittance=HospitalAdmittanceEnum.No,
)


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
            country = "Bahrain"
            country_exposure = sample.travel_exposure
            if not country_exposure:
                country_exposure = country
            res.append(
                NextstrainSampleMetadataInput(
                    strain=sample.lab_id,
                    date=sample.date_collected.now.strftime("%Y-%m-%d") if sample.date_collected else None,
                    country=country,
                    division=division,
                    area=sample.area_name,
                    country_exposure=country_exposure,
                    division_exposure=division,
                    length=sample.genome_length,
                    age=sample.age,
                    sex=sample.gender,
                    originating_lab="Bahrain Public Health Laboratory",
                    submitting_lab="Bahrain Public Health Laboratory",
                    pangolin_lineage=sample.pangolin_lineage,
                    nationality=sample.nationality,
                    ct_value=sample.ct_value,
                    symptoms=sample.symptoms,
                    comorbitities=",".join([x.description for x in sample.comorbidities]),
                    travel_exposure=sample.travel_exposure,
                    hospital_admittance=sample.hospital_admittance,
                )
            )
    res.append(root_genome)
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
