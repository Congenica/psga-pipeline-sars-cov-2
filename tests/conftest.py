# pylint: disable=redefined-outer-name
import os
import sys
from pathlib import Path
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


# path fixtures
@pytest.fixture
def test_data_path() -> Path:
    return Path(__file__).parent / "test_data"


@pytest.fixture
def primer_schemes_dir() -> str:
    return "primer_schemes"


@pytest.fixture
def samples_dir() -> str:
    return "samples"


@pytest.fixture
def check_metadata_data_path(test_data_path: Path) -> Path:
    return test_data_path / "check_metadata"


@pytest.fixture
def concat_csv_data_path(test_data_path: Path) -> Path:
    return test_data_path / "concat_csv"


@pytest.fixture
def contamination_removal_data_path(test_data_path: Path) -> Path:
    return test_data_path / "contamination_removal"


@pytest.fixture
def fetch_primers_data_path(test_data_path: Path) -> Path:
    return test_data_path / "fetch_primers"


@pytest.fixture
def fetch_primers_primer_schemes_data_path(fetch_primers_data_path: Path, primer_schemes_dir: str) -> Path:
    return fetch_primers_data_path / primer_schemes_dir


@pytest.fixture
def integration_test_validation_data_path(test_data_path: Path) -> Path:
    return test_data_path / "integration_test_validation"


@pytest.fixture
def pipeline_results_files_data_path(test_data_path: Path) -> Path:
    return test_data_path / "pipeline_results_files"


@pytest.fixture
def primer_autodetection_data_path(test_data_path: Path) -> Path:
    return test_data_path / "primer_autodetection"


@pytest.fixture
def primer_autodetection_primer_schemes_data_path(
    primer_autodetection_data_path: Path, primer_schemes_dir: str
) -> Path:
    return primer_autodetection_data_path / primer_schemes_dir


@pytest.fixture
def primer_autodetection_sample_dir_data_path(primer_autodetection_data_path: Path, samples_dir: str) -> Path:
    return primer_autodetection_data_path / samples_dir


@pytest.fixture
def typing_data_path(test_data_path: Path) -> Path:
    return test_data_path / "typing"


@pytest.fixture
def variant_definitions_data_path(test_data_path: Path) -> Path:
    return test_data_path / "variant_definitions"
