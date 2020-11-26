from click.testing import CliRunner


def test_reheader_fasta(tmp_path, qc_fasta_files_dir, scripts_test):
    from reheader_fasta import reheader_fasta, SEQUENCE_DESCRIPTION

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
