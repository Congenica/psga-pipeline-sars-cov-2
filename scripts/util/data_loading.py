from pathlib import Path
from typing import Dict
import json
import yaml


def load_json(json_file: str) -> Dict:
    """
    Load a json file
    :param json_file: the json file to load
    :return: the json file as a dictionary
    """
    with open(json_file) as json_handle:
        return json.load(json_handle)


def load_yaml(input_yaml: Path) -> Dict:
    """
    Load a yaml file to dictionary
    """
    with open(input_yaml, "r", encoding="utf8") as yaml_stream:
        return yaml.safe_load(yaml_stream)


def write_yaml(data: Dict, output_yaml: Path) -> None:
    """
    Write a dictionary to a yaml file
    """
    with open(output_yaml, "w", encoding="utf8") as yaml_stream:
        yaml.dump(data, yaml_stream, explicit_start=True, default_flow_style=False, allow_unicode=True)
