import csv
import re
from pathlib import Path
from typing import List
import click


def get_num_from_str(text: str) -> float:
    """
    Extract a number (int or float) from a string
    """
    num = re.search(r"[-+]?(?:\d*\.\d+|\d+)", text).group()
    return float(num)


def convert_time_to_seconds(custom_time: str) -> float:
    """
    Convert Nextflow time to seconds
    """
    time_split = custom_time.split()
    seconds = 0.0
    for t in time_split:
        if "ms" in t:
            seconds += 0.001 * get_num_from_str(t)
        elif "s" in t:
            seconds += get_num_from_str(t)
        elif "m" in t:
            seconds += 60 * get_num_from_str(t)
        elif "h" in t:
            seconds += 60 * 60 * get_num_from_str(t)
        elif "d" in t:
            seconds += 24 * 60 * 60 * get_num_from_str(t)
        elif "w" in t:
            seconds += 7 * 24 * 60 * 60 * get_num_from_str(t)
        elif t == "0":
            pass
        else:
            raise Exception(f"Not sure what this is: {t}")
    return seconds


def convert_size_to_kb(custom_size: str):
    """
    Convert Nextflow size to kilobytes
    """
    kb = 0.0  # kilobytes
    if " B" in custom_size:
        kb += 0.001 * get_num_from_str(custom_size)
    elif " KB" in custom_size:
        kb += get_num_from_str(custom_size)
    elif " MB" in custom_size:
        kb += 1000 * get_num_from_str(custom_size)
    elif " GB" in custom_size:
        kb += 1000000 * get_num_from_str(custom_size)
    elif " TB" in custom_size:
        kb += 1000000000 * get_num_from_str(custom_size)
    elif custom_size == "0":
        pass
    else:
        raise Exception(f"Not sure what this is: {custom_size}")
    return kb


def convert(input_path: Path, output_path: Path, fieldnames: List[str]) -> None:
    """
    Extract specific columns from the trace input report and convert them to seconds and kilobytes.
    Finally, print a summary.
    """

    total_duration = 0.0
    total_realtime = 0.0
    total_peak_rss = 0.0
    total_peak_vmem = 0.0
    total_rchar = 0.0
    total_wchar = 0.0

    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            converted_dict = {
                "name": row["name"],
                "status": row["status"],
                "duration": convert_time_to_seconds(row["duration"]),
                "realtime": convert_time_to_seconds(row["realtime"]),
                "peak_rss": convert_size_to_kb(row["peak_rss"]),
                "peak_vmem": convert_size_to_kb(row["peak_vmem"]),
                "rchar": convert_size_to_kb(row["rchar"]),
                "wchar": convert_size_to_kb(row["wchar"]),
            }
            total_duration += converted_dict["duration"]
            total_realtime += converted_dict["realtime"]
            total_peak_rss += converted_dict["peak_rss"]
            total_peak_vmem += converted_dict["peak_vmem"]
            total_rchar += converted_dict["rchar"]
            total_wchar += converted_dict["wchar"]

            writer.writerow(converted_dict)

    print("Summary")
    print("duration(s)\trealtime(s)\tpeak_rss(KB)\tpeak_vmem(KB)\trchar(KB)\twchar(KB)")
    print(f"{total_duration}\t{total_realtime}\t{total_peak_rss}\t{total_peak_vmem}\t{total_rchar}\t{total_wchar}")


@click.command()
@click.option(
    "--input-tsv",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="Input TSV file. This is a Nextflow trace report",
)
@click.option(
    "--output-tsv",
    type=click.Path(file_okay=True, writable=True),
    required=True,
    help="The simplified TSV report",
)
def process_nextflow_trace_report(input_tsv: str, output_tsv: str) -> None:
    """
    Receive a nextflow trace report as input file and outputs a simplified report
    with same units per column (s, KB)
    """
    input_path = Path(input_tsv)
    output_path = Path(output_tsv)
    fieldnames = ["name", "status", "duration", "realtime", "peak_rss", "peak_vmem", "rchar", "wchar"]

    convert(input_path, output_path, fieldnames)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    process_nextflow_trace_report()
