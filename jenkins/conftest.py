# pylint: disable=redefined-outer-name
import pytest


def pytest_addoption(parser):
    parser.addoption("--results-csv", action="store", required=True)
    parser.addoption("--expected-results-csv", action="store", required=True)
    parser.addoption("--output-path", action="store", required=True)
    parser.addoption("--pathogen", action="store", required=True)
    parser.addoption("--sequencing-technology", action="store", required=True)


@pytest.fixture(scope="session")
def results_csv(pytestconfig):
    return pytestconfig.getoption("results_csv")


@pytest.fixture(scope="session")
def expected_results_csv(pytestconfig):
    return pytestconfig.getoption("expected_results_csv")


@pytest.fixture(scope="session")
def output_path(pytestconfig):
    return pytestconfig.getoption("output_path")


@pytest.fixture(scope="session")
def pathogen(pytestconfig):
    return pytestconfig.getoption("pathogen")


@pytest.fixture(scope="session")
def sequencing_technology(pytestconfig):
    return pytestconfig.getoption("sequencing_technology")
