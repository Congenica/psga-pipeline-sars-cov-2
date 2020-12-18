#!/usr/bin/env python

import csv

import click
from sqlalchemy import and_

from db.database import session_handler
from db.models import Sample


MICROREACT_FIELDS = [
    "id",
    "latitude",
    "longitude",
    "year",
    "month",
    "day",
    "clade",
    "clade__color",
    "lineage",
    "location",
    "location__color",
    "outcome",
    "outcome__color",
    "age",
    "sex",
]


@click.command()
@click.option(
    "--output", required=True, type=click.File("x"), help="The TSV file to write the Microreact input data to"
)
def generate_microreact_input(output):

    with session_handler() as session:

        writer = csv.DictWriter(output, fieldnames=MICROREACT_FIELDS, delimiter="\t")
        writer.writeheader()

        samples = (
            session.query(Sample)
            .filter(and_(Sample.metadata_loaded, Sample.area))
            .order_by(Sample.date_collected)
            .all()
        )
        for sample in samples:
            tsv_row = {
                "id": sample.lab_id,
                "latitude": sample.area.latitude,
                "longitude": sample.area.longitude,
                "year": sample.date_collected.year,
                "month": sample.date_collected.month,
                "day": sample.date_collected.day,
                "clade": None,
                "clade__color": None,
                "lineage": sample.pangolin_lineage,
                "location": sample.area.name,
                "location__color": "#" + sample.area.colour.hex(),
                "outcome": None,
                "outcome__color": None,
                "age": sample.age,
                "sex": None,
            }
            writer.writerow(tsv_row)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    generate_microreact_input()
