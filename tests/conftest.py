# pylint: disable=redefined-outer-name
import os
import sys
from pathlib import Path
import requests
import pytest
from pytest_socket import disable_socket


def pytest_runtest_setup():
    disable_socket()


@pytest.fixture
def fasta_file_generator():
    # generate consensus fasta files
    def generator(path, extension, content):
        for sample_id in ["1", "ERR12313", "ERR12313_barcode001", "529d82ab-51b7-433c-8be6-f1d59759b943"]:
            sample_id = f"NNN{str(sample_id).zfill(5)}"
            fasta_file = path / f"{sample_id}.consensus.{extension}"
            fasta_file.write_text(content.format(sample_id=sample_id))

    return generator


@pytest.fixture(autouse=True)
def scripts_test():
    scripts_path = os.path.dirname(__file__).split("/")
    scripts_path[-1] = "scripts"
    sys.path.append("/".join(scripts_path))


@pytest.fixture
def test_data_path():
    return Path(__file__).parent / "test_data"


@pytest.fixture
def test_data_path_genbank_input(test_data_path):
    return test_data_path / "genbank_input"


@pytest.fixture
def test_data_path_genbank_reference(test_data_path):
    return test_data_path / "genbank_reference"


@pytest.fixture
def setup_splunk(monkeypatch):
    """
    Mock splunk env vars.
    """
    monkeypatch.setenv("HEC_ENABLED_CHILD", "true")
    monkeypatch.setenv("HOSTNAME", "hostname")
    monkeypatch.setenv("HEC_ENDPOINT", "http://hec:1234")
    monkeypatch.setenv("HEC_TOKEN", "token")
    monkeypatch.setenv("HEC_SOURCE", "hec_source")
    monkeypatch.setenv("HEC_SOURCETYPE", "hec_sourcetype")


# pylint: disable=unused-argument
def mocked_request(*args, **kwargs):
    class MockResponse:
        def __init__(self, json_data, status_code, url):
            self.json_data = json_data
            self.status_code = status_code
            self.url = url

        def json(self):
            return self.json_data

        def raise_for_status(self):
            """I basically copy pasted this from the real raise_for_status"""

            reason = "MOCK RAISE"
            http_error_msg = ""
            if 400 <= self.status_code < 500:
                http_error_msg = f"{self.status_code} Client Error: {reason} for url: {self.url}"

            elif 500 <= self.status_code < 600:
                http_error_msg = f"{self.status_code} Server Error: {reason} for url: {self.url}"

            if http_error_msg:
                raise requests.exceptions.HTTPError(http_error_msg, response=self)

    url_data = {
        "auth": [{}, 200],
        "service_auth": [{}, 200],
    }
    url = args[1]  # request.sessions.Session

    return MockResponse(url_data, 200, url=url)


@pytest.fixture(autouse=True)
def no_requests(monkeypatch):
    """
    Ensure no real requests get sent,
    overriding get/post/put with a mocked version
    """

    # pylint: disable=unused-argument
    def fail(*args, **kwargs):
        raise Exception("An actual requests.sessions.Session.request was called from a test")

    # don't let someone accidentally send a request here
    monkeypatch.setattr("requests.sessions.Session.request", fail)
    monkeypatch.setattr("requests.get", mocked_request)
    monkeypatch.setattr("requests.put", mocked_request)
    monkeypatch.setattr("requests.post", mocked_request)

    monkeypatch.setattr("requests.sessions.Session.get", mocked_request)
    monkeypatch.setattr("requests.sessions.Session.put", mocked_request)
    monkeypatch.setattr("requests.sessions.Session.post", mocked_request)
