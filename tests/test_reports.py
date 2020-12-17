from collections import namedtuple
from datetime import datetime, timedelta

import pytest


@pytest.mark.parametrize(
    "mocked_response",
    [
        lambda content: [False, b""],  # failed fetch, testing local data
        lambda content: [True, content.read_bytes()],  # successful data pull
    ],
)
def test_report_strain_level_and_global_context(
    sample_generator, tmp_path, monkeypatch, test_data_pangolearn, mocked_response, report_generator
):
    def mock_response(session_obj, url):
        content = test_data_pangolearn / url.split("/")[-1]
        response = namedtuple("MockResponse", ["ok", "content"])
        return response(*mocked_response(content))

    monkeypatch.setattr("requests.sessions.Session.get", mock_response)

    sample_generator()

    result, report = report_generator(
        "strain_level_and_global_context",
        "--lineage-notes-url",
        "http://test.url/lineage_notes.txt",
        "--metadata-url",
        "http://test.url/lineages.metadata.csv",
        "--pangolearn-dir",
        str(test_data_pangolearn),
    )

    assert result.exit_code == 0
    assert report == [
        [
            "B.1",
            "A large European lineage that corresponds to the Italian outbreak.",
            "1",
            "100%",
            "2",
        ],
        ["", "", "1", "", "2"],
    ]


def test_report_strain_prevalence(sample_generator, report_generator):
    area_name = "JERDAB"
    governorate_name = "Capital"
    prevalent_strain = "B.1.366"

    for lineage in (prevalent_strain, "B.1.369"):
        sample_generator(governorate_name, area_name, lineage, datetime.now())
        sample_generator(governorate_name, area_name, lineage, datetime.now() - timedelta(weeks=1))
    sample_generator(governorate_name, area_name, prevalent_strain, datetime.now())

    result, report = report_generator("strain_prevalence")

    assert result.exit_code == 0

    totals = [
        "0.67",  # 2/3 of records
        "0.5",  # 1/2 in last week
        "0.17",  # 2/3 - 1/2 difference
        "2",  # prevalent_strain has 2 records in (today - 1 week)
    ]
    assert report[0] == [area_name, prevalent_strain] + totals
    assert report[1] == [governorate_name, prevalent_strain] + totals
    assert report[2] == ["Bahrain", prevalent_strain] + totals


def test_report_strain_first_seen(sample_generator, report_generator):
    date = "2020-12-17"
    sample_generator(date_collected=date)

    result, report = report_generator("strain_first_seen")

    assert result.exit_code == 0
    assert report == [["B.1", date, "", "", "JERDAB", "", "", "0"]]


def test_report_sample_dump(sample_generator, report_generator):
    sample_generator()

    result, report = report_generator("sample_dump")

    assert result.exit_code == 0
    assert len(report) == 1  # single record
