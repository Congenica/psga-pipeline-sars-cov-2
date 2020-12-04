from os.path import join

from click.testing import CliRunner


def test_fasta_concatenation(tmp_path, fasta_file_generator, root_genome):
    from reheader_fasta import SEQUENCE_DESCRIPTION, FASTA_FILE_HANDLE
    from concatenate_fasta import concatenate_fasta

    fasta_file_generator(
        path=tmp_path,
        extension=FASTA_FILE_HANDLE,
        content=">Consensus_{sample_id}" + f"{SEQUENCE_DESCRIPTION}\n{'N' * 29903}\n",
    )
    headers = [fasta_file.read_text().split("\n")[0] for fasta_file in tmp_path.glob("*")]
    result_filename = "concatenated.fasta"

    args = ["--output", join(tmp_path, result_filename), "--root-genome", str(root_genome), str(tmp_path)]

    runner = CliRunner()
    result = runner.invoke(concatenate_fasta, args)

    assert result.exit_code == 0
    concatenated_fasta = next(tmp_path.glob(result_filename)).read_text()
    # check file headers in concatenated file
    assert all(header in concatenated_fasta for header in headers)
    # check SARS-CoV-2 prepended
    assert concatenated_fasta.startswith(f">NC_045512.2 {SEQUENCE_DESCRIPTION}")
