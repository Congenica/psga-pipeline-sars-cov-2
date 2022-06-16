import pytest
from click.testing import CliRunner

from scripts.sars_cov_2.reheader_fasta import reheader_fasta, FASTA_FILE_EXTENSION, FASTA_FILE_HANDLE


@pytest.mark.parametrize(
    "input_fasta_header",
    [
        ">Consensus_{sample_id}",
        ">{sample_id}",
        ">blablabla{sample_id}",
        ">{sample_id}blablabla",
        ">{sample_id}.blablabla",
        ">{sample_id}_blablabla",
        ">{sample_id} - Calmodulin - Human, rabbit, bovine, rat, and chicken",
        ">{sample_id}|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]",
    ],
)
def test_reheader_fasta(tmp_path, fasta_file_generator, input_fasta_header):

    fasta_file_generator(
        path=tmp_path,
        extension=FASTA_FILE_EXTENSION,
        content=f"{input_fasta_header}\n{'N' * 29903}\n",
    )

    runner = CliRunner()
    result = runner.invoke(
        reheader_fasta,
        ["--input-dir", tmp_path, "--output-dir", tmp_path],
    )

    assert result.exit_code == 0

    result_files = list(tmp_path.glob(f"*.{FASTA_FILE_HANDLE}"))
    assert len(result_files) == len(list(tmp_path.glob(f"*.{FASTA_FILE_EXTENSION}")))

    for reheadered_fasta in result_files:
        sample_id = reheadered_fasta.name.replace(f".{FASTA_FILE_HANDLE}", "")
        content = reheadered_fasta.read_text()
        assert content.startswith(f">{sample_id}\n")