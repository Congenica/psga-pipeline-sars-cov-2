from pathlib import Path
from click.testing import CliRunner

from scripts.sars_cov_2.extract_typing_data import extract_typing_data
from scripts.util.data_loading import load_yaml


def test_extract_typing_data(tmp_path: Path, typing_data_path: Path, variant_definitions_data_path: Path):

    output_yaml_path = tmp_path / "definitions.yaml"
    expected_output_yaml_path = typing_data_path / "expected_sars_cov_2_types.yaml"

    rv = CliRunner().invoke(
        extract_typing_data,
        [
            "--yaml-input-dir",
            variant_definitions_data_path,
            "--output-yaml-file",
            output_yaml_path,
        ],
    )
    assert rv.exit_code == 0
    assert load_yaml(output_yaml_path) == load_yaml(expected_output_yaml_path)
