#!/usr/bin/env python

import csv
from datetime import date
import re

import click
from click import ClickException
from sqlalchemy import String, cast, func

from db.database import session_handler
from db.models import Area, Governorate, Sample

EXPECTED_HEADERS = [
    "MRN",
    "AGE",
    "NATIONALITY",
    "GOVERNORATE",
    "AREA",
    "BLOCK",
    "SAMPLE ID",
    "ASSIGN DATE",
    "CT",
    "SYMPTOMS",
]


def _validate_and_normalise_row(session, row):
    # strip leading and trailing spaces from everything
    for f in EXPECTED_HEADERS:
        row[f] = row[f].lstrip().rstrip()

    errs = []

    # MRN should be only digits
    if re.match(r"\d+$", row["MRN"]):
        row["MRN"] = int(row["MRN"])
    else:
        errs.append(f"MRN \"{row['MRN']}\" is not an integer")

    # AGE should be digits, space, "Y". turn into just digits. if not Y make it 0
    if age := re.match(r"(\d{1,3}) ([A-Za-z])$", row["AGE"]):
        if age.group(2).lower() == "y":
            row["AGE"] = int(age.group(1))
        else:
            row["AGE"] = 0
    else:
        errs.append(f"AGE \"{row['AGE']}\" should be a number followed by a space then a letter (probably Y)")

    # NATIONALITY should be str
    if not re.match(r"[\w ]+$", row["NATIONALITY"]):
        errs.append(f"NATIONALITY \"{row['NATIONALITY']}\" doesn't look like words")

    # GOVERNORATE should be in governorate enum
    if governorate := re.match(r"([\w ]+)$", row["GOVERNORATE"]):
        found_governorate = (
            session.query(Governorate)
            .filter(func.lower(cast(Governorate.name, String)) == governorate.group(1).lower())
            .one_or_none()
        )
        if found_governorate:
            row["GOVERNORATE"] = found_governorate
        else:
            errs.append(f"GOVERNORATE \"{row['GOVERNORATE']}\" looked right, but not found in database")
    else:
        errs.append(f"GOVERNORATE \"{row['GOVERNORATE']}\" doesn't look like words")

    # AREA should be in the area table
    if area := re.match(r"([\w /']+)$", row["AREA"]):
        # do a case insensitive match just in case, also the spaces around punctuation is inconsistent, so we compare
        # ignoring whitespace
        whitespace = re.compile(r"\s+")
        normalised_area = whitespace.sub("", area.group(1).lower())
        found_area = (
            session.query(Area)
            .filter(func.regexp_replace(func.lower(Area.name), r"\s+", "", "g") == normalised_area)
            .one_or_none()
        )
        if found_area:
            row["AREA"] = found_area
        else:
            errs.append(f"AREA \"{row['AREA']}\" looked right, but not found in database")
    else:
        errs.append(f"AREA \"{row['AREA']}\" doesn't look like words")

    # BLOCK should be a number, appears to normally be 3 digits, but have also seen 4, so we'll be lax on this
    if not re.match(r"(\d+)$", row["BLOCK"]):
        errs.append(f"BLOCK \"{row['BLOCK']}\" should contain digits")

    # SAMPLE ID should be numbers, but not an integer
    if not re.match(r"\d+$", row["SAMPLE ID"]):
        errs.append(f"SAMPLE ID \"{row['SAMPLE ID']}\" is not a number")

    # ASSIGN DATE should be dd/mm/yyyy
    if assign_date_match := re.match(r"(\d{2})/(\d{2})/(\d{4})$", row["ASSIGN DATE"]):
        try:
            row["ASSIGN DATE"] = date(
                int(assign_date_match.group(3)), int(assign_date_match.group(2)), int(assign_date_match.group(1))
            )
        except ValueError as e:
            errs.append(f"ASSIGN DATE \"{row['ASSIGN DATE']}\" is not a valid date: {e}")
    else:
        errs.append(f"ASSIGN DATE \"{row['ASSIGN DATE']}\" doesn't look like a date")

    # CT should be CT: followed by a float, which might not have the fractional part
    # it's also not in every row in the example data, so assuming it's optional
    if ct := re.match(r"CT:(\d+(?:\.\d+)?)$", row["CT"]):
        row["CT"] = float(ct.group(1))
    else:
        row["CT"] = None

    # SYMPTOMS doesn't really need validating

    if errs:
        raise ValueError("\n".join(errs))

    return row


@click.command()
@click.option("--file", required=True, type=click.File("r"), help="The metadata TSV input file")
def load_iseha_data(file):
    """
    Read in a TSV file of I-SEHA metadata and load each row into the database. Invalid rows are warned about, but
    skipped over.

    :param File file: the TSV file to load
    """
    reader = csv.DictReader(file, delimiter="\t")

    found_headers = reader.fieldnames
    for expected_header in EXPECTED_HEADERS:
        if expected_header not in found_headers:
            err = (
                "Could not find expected header '"
                + expected_header
                + "'. Got:\n"
                + ", ".join(found_headers)
                + "\nexpected\n"
                + ", ".join(EXPECTED_HEADERS)
            )
            raise ClickException(err)

    inserted = set()
    updated = set()
    errors = set()
    with session_handler() as session:
        for row in reader:
            try:
                row = _validate_and_normalise_row(session, row)
                if existing_sample := session.query(Sample).filter(Sample.lab_id == row["SAMPLE ID"]).one_or_none():
                    existing_sample.mrn = row["MRN"]
                    existing_sample.age = row["AGE"]
                    existing_sample.nationality = row["NATIONALITY"]
                    existing_sample.governorate = row["GOVERNORATE"]
                    existing_sample.area = row["AREA"]
                    existing_sample.block_number = row["BLOCK"]
                    existing_sample.date_collected = row["ASSIGN DATE"]
                    existing_sample.ct_value = row["CT"]
                    existing_sample.symptoms = row["SYMPTOMS"]
                    updated.add(row["SAMPLE ID"])
                else:
                    sample = Sample(
                        mrn=row["MRN"],
                        age=row["AGE"],
                        nationality=row["NATIONALITY"],
                        governorate=row["GOVERNORATE"],
                        area=row["AREA"],
                        block_number=row["BLOCK"],
                        lab_id=row["SAMPLE ID"],
                        sample_number=int(row["SAMPLE ID"]),
                        date_collected=row["ASSIGN DATE"],
                        ct_value=row["CT"],
                        symptoms=row["SYMPTOMS"],
                    )
                    session.add(sample)
                    inserted.add(row["SAMPLE ID"])
                session.commit()
            except ValueError as e:
                click.echo(f"Invalid row for sample ID {row['SAMPLE ID']}:\n{e}", err=True)
                errors.add(row["SAMPLE ID"])

    click.echo("Inserted samples: " + ", ".join(map(str, inserted)))
    click.echo("Updated samples: " + ", ".join(map(str, updated)))
    if errors:
        click.echo(
            "Errors encountered, these samples were not inserted or updated: " + ", ".join(map(str, errors)), err=True
        )


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_iseha_data()
