# pylint: disable=redefined-outer-name
import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from reheader_fasta import reheader_fasta, SEQUENCE_DESCRIPTION


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


def test_reheader_fasta(tmp_path, qc_fasta_files_dir):
    args = [str(qc_fasta_files_dir), str(tmp_path)]
    runner = CliRunner()
    result = runner.invoke(reheader_fasta, args)

    assert result.exit_code == 0

    result_files = list(tmp_path.glob("*.fasta"))
    assert len(result_files) == len(list(qc_fasta_files_dir.glob("*.fa")))

    for reheadered_fasta in result_files:
        sample_id = reheadered_fasta.name.replace(".fasta", "")
        content = reheadered_fasta.read_text()
        assert content.startswith(f">{sample_id} {SEQUENCE_DESCRIPTION}\n")
