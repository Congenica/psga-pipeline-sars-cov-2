from click.testing import CliRunner

from scripts.sars_cov_2.extract_typing_data import extract_typing_data
from scripts.util.data_loading import load_yaml


def test_extract_typing_data(tmp_path, test_data_path):

    input_path = test_data_path / "variant_definitions"
    output_yaml_path = tmp_path / "definitions.yaml"
    expected_output_yaml_path = test_data_path / "expected_sars_cov_2_types.yaml"

    rv = CliRunner().invoke(
        extract_typing_data,
        [
            "--yaml-input-dir",
            input_path,
            "--output-yaml-file",
            output_yaml_path,
        ],
    )
    assert rv.exit_code == 0
    assert load_yaml(output_yaml_path) == load_yaml(expected_output_yaml_path)
