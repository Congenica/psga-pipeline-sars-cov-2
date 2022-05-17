from typing import Dict
import json


def load_json(json_file: str) -> Dict:
    """
    Load a json file
    :param json_file: the json file to load
    :return: the json file as a dictionary
    """
    with open(json_file) as json_handle:
        json_dict = json.load(json_handle)
    return json_dict
