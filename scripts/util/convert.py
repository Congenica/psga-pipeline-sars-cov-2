from pathlib import Path
import csv
import ast
import json


def csv_to_json(csv_path: Path, json_path: Path, sample_id: str):
    """
    Convert a CSV file to a JSON file using sample_id as key in the JSON
    """
    data = {}

    with open(csv_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if sample_id not in row:
                raise KeyError(f"key {sample_id} was not found in row {row}")
            sample_id_key = row[sample_id]

            # covert strings representing a dictionary to Python dictionaries
            record = {}
            for key, val in row.items():
                if not key:
                    raise ValueError(f"Found a None key in CSV: {csv_path}. Is the CSV malformed?")

                if key == sample_id:
                    # discard as sample_id is the main key of this record
                    continue

                if val and val.startswith("{") and val.endswith("}"):
                    record[key] = ast.literal_eval(val)
                else:
                    record[key] = val

            data[sample_id_key] = record

    with open(json_path, "w", encoding="utf-8") as jsonfile:
        jsonfile.write(json.dumps(data, indent=4, sort_keys=True))
