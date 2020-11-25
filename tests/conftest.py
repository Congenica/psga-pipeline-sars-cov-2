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
def qc_fasta_files_dir(tmp_path):
    for sample_id in range(5):
        sample_id = f"NNN{str(sample_id).zfill(5)}"
        fasta_file = tmp_path / f"{sample_id}.fa"
        fasta_file.write_text(
            f">Consensus_{sample_id}.primertrimmed.consensus_threshold_0.75_quality_20\n" f"{'N' * 29903}\n"
        )
    return tmp_path


@pytest.fixture
def scripts_test():
    scripts_path = os.path.dirname(__file__).split("/")
    scripts_path[-1] = "scripts"
    sys.path.append("/".join(scripts_path))
