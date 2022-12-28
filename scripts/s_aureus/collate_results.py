import os
import json
import csv
import click


@click.command()
@click.option(
    "--metadata-file",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
    help="the sample metadata file",
)
def collate_results(metadata_file):
    # Merge the resultfiles.json files
    output_resultfiles_json_dict = dict()
    for f in os.listdir():
        if "resultfiles.json" not in f or str(f) == "resultfiles.json":
            print("Skipping: " + str(f))
            continue
        print("Loading: " + str(f))
        with open(f) as channel_json_file:
            channel_dict = json.load(channel_json_file)
            output_resultfiles_json_dict.update(channel_dict)

    # Find duplicate files (Framework will error if 2 runs point to the same file)
    all_files = set()
    dup_files = set()
    for sam_key, sam_body in output_resultfiles_json_dict.items():
        for i in sam_body:
            if i["file"] in all_files:
                dup_files.add(i["file"])
            all_files.add(i["file"])
    if dup_files:
        print("WARNING found duplicate files will remove them: " + str(dup_files))

    # Strip the duplicate files out
    new_output_resultfiles_json_dict = dict()
    for sam_key, sam_body in output_resultfiles_json_dict.items():
        new_sam_body = list()
        for i in sam_body:
            if i["file"] not in list(dup_files):
                new_sam_body.append(i)
        new_output_resultfiles_json_dict[sam_key] = new_sam_body

    with open("resultfiles.json", "w") as output_result_files:
        json.dump(new_output_resultfiles_json_dict, output_result_files)

    # Merge the results.csv files
    output_results_csv_list = list()
    for f in os.listdir():
        if "results.csv" not in f or str(f) == "results.csv":
            print("Skipping: " + str(f))
            continue
        print("Loading: " + str(f))
        with open(f) as channel_csv_file:
            # All of these channel results csv files are single sample
            output_results_csv_list.append(next(csv.DictReader(channel_csv_file)))

    # Add on any failed samples where the container failed completely and did not produce a results.csv for that sample.
    # A sample is failed if there are no results in output_results_csv_list, but sample is in the input metadata file
    samples_with_output = set(r["SAMPLE_ID"] for r in output_results_csv_list)
    with open(metadata_file) as input_metadata_file:
        for md_sample in csv.DictReader(input_metadata_file):
            if md_sample["SAMPLE_ID"] not in samples_with_output:
                output_results_csv_list.append({"SAMPLE_ID": md_sample["SAMPLE_ID"], "STATUS": "Failed"})

    with open("results.csv", "w", encoding="utf8") as output_csv_file:
        # Need all keys for all fields in all rows preserve order
        all_fields_ordered = list()
        for row in output_results_csv_list:
            for c in row.keys():
                if c not in all_fields_ordered:
                    all_fields_ordered.append(c)
        dict_writer = csv.DictWriter(output_csv_file, fieldnames=all_fields_ordered)
        dict_writer.writeheader()
        dict_writer.writerows(output_results_csv_list)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    collate_results()
