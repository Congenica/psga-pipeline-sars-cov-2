import csv
from pathlib import Path
from typing import List


def write_list_to_file(elements: List[str], output_path: Path) -> None:
    """
    Store a list of strings to CSV file
    """
    with open(output_path, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=",")
        for element in elements:
            writer.writerow([element])
