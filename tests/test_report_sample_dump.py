import csv

from click.testing import CliRunner
from scripts.generate_report import generate_report


def test_report_sample_dump(sample_generator, tmp_path):
    sample_generator()

    report = tmp_path / "sample_dump.csv"
    result = CliRunner().invoke(generate_report, ["--output", str(report), "--report", "sample_dump"])

    assert result.exit_code == 0

    with open(report, mode="r") as report:
        report = [row for row in csv.reader(report)]

    # header and single record
    assert len(report) == 2
