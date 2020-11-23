from click.testing import CliRunner


def test_fasta_concatenation(tmp_path, fasta_file_generator, scripts_test):
    from reheader_fasta import SEQUENCE_DESCRIPTION, FASTQ_FILE_HANDLE
    from concatenate_fasta import concatenate_fasta, CONCATENATED_FASTA_FILENAME

    fasta_file_generator(
        path=tmp_path,
        extension=FASTQ_FILE_HANDLE,
        content=">Consensus_{sample_id}" + f"{SEQUENCE_DESCRIPTION}\n{'N' * 29903}\n",
    )
    headers = [fasta_file.read_text().split("\n")[0] for fasta_file in tmp_path.glob("*")]

    runner = CliRunner()
    result = runner.invoke(concatenate_fasta, [str(tmp_path)])

    concatenated_fasta = next(tmp_path.glob(CONCATENATED_FASTA_FILENAME)).read_text()
    assert result.exit_code == 0
    # check file headers in concatenated file
    assert all(header in concatenated_fasta for header in headers)
    # check SARS-CoV-2 prepended
    assert concatenated_fasta.startswith(f">Wuhan-Hu-1 {SEQUENCE_DESCRIPTION}")
