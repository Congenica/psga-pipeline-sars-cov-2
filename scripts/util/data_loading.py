from pathlib import Path, PosixPath
from typing import Dict
import json
from json import JSONEncoder
import yaml


class PathJSONEncoder(JSONEncoder):
    """
    Enable JSON serialisation of PosixPath objects
    """

    def default(self, o):
        if isinstance(o, PosixPath):
            return str(o)
        return super().default(o)


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


def write_json(data: Dict, output_json: Path) -> None:
    """
    Write a dictionary to a json file
    """
    with open(output_json, "w", encoding="utf8") as json_stream:
        json.dump(data, json_stream, cls=PathJSONEncoder, sort_keys=True, indent=4)


def write_yaml(data: Dict, output_yaml: Path) -> None:
    """
    Write a dictionary to a yaml file
    """
    with open(output_yaml, "w", encoding="utf8") as yaml_stream:
        yaml.dump(data, yaml_stream, explicit_start=True, default_flow_style=False, allow_unicode=True)
