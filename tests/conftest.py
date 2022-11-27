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


@pytest.fixture
def test_data_path():
    return Path(__file__).parent / "test_data"
