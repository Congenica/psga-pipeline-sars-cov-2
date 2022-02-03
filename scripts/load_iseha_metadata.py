#!/usr/bin/env python

from typing import List
import csv
from datetime import date
import re
from pathlib import Path

import click
from click import ClickException
from sqlalchemy import String, cast, func

from scripts.db.database import session_handler
from scripts.db.models import Area, Governorate, Sample, SampleQC

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
    def location_error_message(location_type, location):
        return (
            f'{location_type} "{location}" is not a recognised location - '
            "please check the list of known locations and for typo's etc"
        )

    # strip leading and trailing spaces from everything
    for f in EXPECTED_HEADERS:
        row[f] = row[f].lstrip().rstrip() if row[f] is not None else ""

    errs = []

    # MRN contains digits. If the person is a foreigner, a 'T' is prefixed
    if not re.match(r"(T)?\d+$", row["MRN"]):
        errs.append(f"MRN \"{row['MRN']}\" is expected to be an integer, prefixed by 'T' if the person is a foreigner")

    # AGE can come in two different formats:
    # AGE can be digits, space, "Y". turn into just digits. if not Y make it 0
    # AGE can be digits-only
    age_column_value = row["AGE"]
    age_y = re.match(r"(\d{1,3}) ([A-Za-z])$", age_column_value)
    if age_y:
        if age_y.group(2).lower() == "y":
            row["AGE"] = int(age_y.group(1))
        else:
            row["AGE"] = 0
    else:
        age_digit = re.match(r"(\d{1,3})$", age_column_value)
        if age_digit:
            row["AGE"] = int(age_digit.group(1))
        else:
            errs.append(
                f"AGE \"{row['AGE']}\" does not follow any known age formats. It needs to be a number-only "
                f"or a number followed by a space then a letter (probably Y)"
            )

    # NATIONALITY should be str
    if not re.match(r"[\w ]+$", row["NATIONALITY"]):
        errs.append(f"NATIONALITY \"{row['NATIONALITY']}\" doesn't look like words")

    # GOVERNORATE should be in governorate enum
    governorate = re.match(r"([\w ]+)$", row["GOVERNORATE"])
    if governorate:
        found_governorate = (
            session.query(Governorate)
            .filter(func.lower(cast(Governorate.name, String)) == governorate.group(1).lower())
            .one_or_none()
        )
        if found_governorate:
            row["GOVERNORATE"] = found_governorate
        else:
            errs.append(location_error_message("GOVERNORATE", row["GOVERNORATE"]))
    elif row["GOVERNORATE"]:
        errs.append(location_error_message("GOVERNORATE", row["GOVERNORATE"]))

    # AREA should be in the area table
    area = re.match(r"([\w /']+)$", row["AREA"])
    if area:
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
            errs.append(location_error_message("AREA", row["AREA"]))
    elif row["AREA"]:
        errs.append(location_error_message("AREA", row["AREA"]))

    # BLOCK should be a number, appears to normally be 3 digits, but have also seen 4, so we'll be lax on this
    if row["BLOCK"] and not re.match(r"(\d+)$", row["BLOCK"]):
        errs.append(f"BLOCK \"{row['BLOCK']}\" should contain digits")

    # ASSIGN DATE should be dd/mm/yyyy
    assign_date_match = re.match(r"(\d{2})/(\d{2})/(\d{4})$", row["ASSIGN DATE"])
    if assign_date_match:
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
    ct = re.match(r"CT:(\d+(?:\.\d+)?)$", row["CT"])
    if ct:
        row["CT"] = float(ct.group(1))
    else:
        row["CT"] = None

    # SYMPTOMS doesn't really need validating

    if errs:
        raise ValueError("\n".join(errs))

    return row


def write_list_to_file(out_file: Path, elements: List[str]) -> None:
    with open(out_file, "w") as f:
        for element in elements:
            f.write(f"{element}\n")


def write_sample_list_files(
    updated_samples: List[str],
    valid_samples: List[str],
    file_all_samples: str,
    file_current_samples: str,
    file_updated_samples: str,
    file_qc_pass_samples: str,
) -> None:
    if file_all_samples:
        with session_handler() as session:
            all_samples_with_metadata = session.query(Sample.lab_id).filter(Sample.metadata_loaded).all()
            all_samples_with_metadata = [x[0] for x in all_samples_with_metadata]
            write_list_to_file(Path(file_all_samples), all_samples_with_metadata)
    if file_current_samples:
        write_list_to_file(Path(file_current_samples), valid_samples)
    if file_qc_pass_samples:
        with session_handler() as session:
            samples_with_qc_pass = session.query(Sample.lab_id).join(Sample.sample_qc).filter(SampleQC.qc_pass).all()
            samples_with_qc_pass = [x[0] for x in samples_with_qc_pass]
            write_list_to_file(Path(file_qc_pass_samples), samples_with_qc_pass)
    if file_updated_samples:
        write_list_to_file(Path(file_updated_samples), updated_samples)


@click.command()
@click.option("--file", required=True, type=click.File("r"), help="The metadata TSV input file")
@click.option(
    "--output_all_samples_with_metadata",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with all ALL sample names found historically, "
    "which have metadata loaded in the database",
)
@click.option(
    "--output_current_samples_with_metadata",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with sample names within the current metadata file, "
    "which were loaded into the database",
)
@click.option(
    "--output_samples_with_qc_pass",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with all sample names, "
    "which were processed by artic-ncov2019 pipeline and were marked as QC_PASS",
)
@click.option(
    "--output_samples_updated",
    required=False,
    type=click.Path(file_okay=True, writable=True),
    help="File path to populate with all sample names, which were re-submitted and overwritten with provided metadata",
)
def load_iseha_data(
    file,
    output_all_samples_with_metadata,
    output_current_samples_with_metadata,
    output_samples_with_qc_pass,
    output_samples_updated,
):
    """
    Read in a TSV file of I-SEHA metadata and load each row into the database. Invalid rows are warned about, but
    skipped over.
    """
    reader = csv.DictReader(file, delimiter="\t")

    headers = reader.fieldnames
    if set(headers) != set(EXPECTED_HEADERS):
        err = "Unexpected TSV headers, got:\n" + ", ".join(headers) + "\nexpected\n" + ", ".join(EXPECTED_HEADERS)
        raise ClickException(err)

    file_samples_with_valid_meta = []
    file_samples_updated = []

    inserted = set()
    updated = set()
    errors = set()
    with session_handler() as session:
        for row in reader:
            try:
                row = _validate_and_normalise_row(session, row)
                governorate = row["GOVERNORATE"] if row["GOVERNORATE"] else None
                area = row["AREA"] if row["AREA"] else None
                block = row["BLOCK"] if row["BLOCK"] else None
                sample_name = row["SAMPLE ID"]
                existing_sample = session.query(Sample).filter(Sample.lab_id == sample_name).one_or_none()
                if existing_sample:
                    existing_sample.mrn = row["MRN"]
                    existing_sample.age = row["AGE"]
                    existing_sample.nationality = row["NATIONALITY"]
                    existing_sample.governorate = governorate
                    existing_sample.area = area
                    existing_sample.block_number = block
                    existing_sample.date_collected = row["ASSIGN DATE"]
                    existing_sample.ct_value = row["CT"]
                    existing_sample.symptoms = row["SYMPTOMS"]
                    existing_sample.metadata_loaded = True
                    updated.add(sample_name)
                    file_samples_updated.append(sample_name)
                else:
                    sample = Sample(
                        mrn=row["MRN"],
                        age=row["AGE"],
                        nationality=row["NATIONALITY"],
                        lab_id=sample_name,
                        sample_number=sample_name,
                        governorate=governorate,
                        area=area,
                        block_number=block,
                        date_collected=row["ASSIGN DATE"],
                        ct_value=row["CT"],
                        symptoms=row["SYMPTOMS"],
                        metadata_loaded=True,
                    )
                    session.add(sample)
                    inserted.add(sample_name)
                session.commit()
                file_samples_with_valid_meta.append(sample_name)
            except ValueError as e:
                click.echo(f"Invalid row for sample ID {row['SAMPLE ID']}:\n{e}", err=True)
                errors.add(row["SAMPLE ID"])

    write_sample_list_files(
        updated_samples=file_samples_updated,
        valid_samples=file_samples_with_valid_meta,
        file_all_samples=output_all_samples_with_metadata,
        file_current_samples=output_current_samples_with_metadata,
        file_updated_samples=output_samples_updated,
        file_qc_pass_samples=output_samples_with_qc_pass,
    )

    if errors:
        raise ClickException("Errors encountered: " + ", ".join(map(str, errors)))
    click.echo("Inserted samples: " + ", ".join(map(str, inserted)))
    click.echo("Updated samples: " + ", ".join(map(str, updated)))


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    load_iseha_data()
