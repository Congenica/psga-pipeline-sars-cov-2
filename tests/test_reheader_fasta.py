from click.testing import CliRunner


def test_reheader_fasta(tmp_path, fasta_file_generator, scripts_test):
    from reheader_fasta import reheader_fasta, SEQUENCE_DESCRIPTION, FASTQ_FILE_EXTENSION, FASTQ_FILE_HANDLE
    fasta_file_generator(
        path=tmp_path,
        extension=FASTQ_FILE_EXTENSION,
        content=">Consensus_{sample_id}" + f"{SEQUENCE_DESCRIPTION}\n{'N' * 29903}\n",
    )

    runner = CliRunner()
    result = runner.invoke(reheader_fasta, [str(tmp_path)])

    assert result.exit_code == 0

    result_files = list(tmp_path.glob(f"*.{FASTQ_FILE_HANDLE}"))
    assert len(result_files) == len(list(tmp_path.glob(f"*.{FASTQ_FILE_EXTENSION}")))

    for reheadered_fasta in result_files:
        sample_id = reheadered_fasta.name.replace(f".{FASTQ_FILE_HANDLE}", "")
        content = reheadered_fasta.read_text()
        assert content.startswith(f">{sample_id} {SEQUENCE_DESCRIPTION}\n")
