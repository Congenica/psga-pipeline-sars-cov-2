import csv
from datetime import datetime, timedelta

from click.testing import CliRunner
from scripts.generate_report import generate_report


def test_report_strain_prevalence(sample_generator, tmp_path):
    area_name = "JERDAB"
    governorate_name = "Capital"

    prevalent_strain = "B.1.366"
    for lineage in (prevalent_strain, "B.1.369"):
        sample_generator(governorate_name, area_name, lineage, datetime.now())
        sample_generator(governorate_name, area_name, lineage, datetime.now() - timedelta(weeks=1))
    sample_generator(governorate_name, area_name, prevalent_strain, datetime.now())

    runner = CliRunner()
    report = tmp_path / "strain_prevalence.csv"
    result = runner.invoke(generate_report, ["--output", str(report), "--report", "strain_prevalence"])

    assert result.exit_code == 0

    with open(report, mode="r") as report:
        report = [row for row in csv.reader(report)]

    totals = [
        "0.67",  # 2/3 of records
        "0.5",  # 1/2 in last week
        "0.17",  # 2/3 - 1/2 difference
        "2",  # prevalent_strain has 2 records in (today - 1 week)
    ]

    assert report[1] == [area_name, prevalent_strain] + totals
    assert report[2] == [governorate_name, prevalent_strain] + totals
    assert report[3] == ["Bahrain", prevalent_strain] + totals
