# pylint: disable=redefined-outer-name
import os
import sys
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def tmp_path():
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def fasta_file_generator():
    def generator(path, extension, content):
        for sample_id in range(5):
            sample_id = f"NNN{str(sample_id).zfill(5)}"
            fasta_file = path / f"{sample_id}.{extension}"
            fasta_file.write_text(content.format(sample_id=sample_id))

    return generator


@pytest.fixture(autouse=True)
def scripts_test():
    scripts_path = os.path.dirname(__file__).split("/")
    scripts_path[-1] = "scripts"
    sys.path.append("/".join(scripts_path))


@pytest.fixture
def root_genome():
    return Path(__file__).parent.parent / "data" / "FASTA" / "SARS-CoV-2.fasta"
