import os
import json
import csv

def collate_results():
    # Merge the resultfiles.json files
    output_resultfiles_json_dict = dict()
    for f in os.listdir():
        if 'resultfiles.json' not in f or str(f) == 'resultfiles.json':
            print('Skipping: ' + str(f))
            continue
        print('Loading: ' + str(f))
        with open(f) as channel_json_file:
            channel_dict = json.load(channel_json_file)
            output_resultfiles_json_dict.update(channel_dict)

    with open('resultfiles.json', 'w') as output_result_files:
        json.dump(output_resultfiles_json_dict, output_result_files)

    # Merge the results.csv files
    output_results_csv_list = list()
    for f in os.listdir():
        if 'results.csv' not in f or str(f) == 'results.csv':
            print('Skipping: ' + str(f))
            continue
        print('Loading: ' + str(f))
        with open(f) as channel_csv_file:
            # All of these channel results csv files are single sample
            output_results_csv_list.append(next(csv.DictReader(channel_csv_file)))

    with open('results.csv', 'w', encoding='utf8') as output_csv_file:
        dict_writer = csv.DictWriter(output_csv_file, fieldnames=output_results_csv_list[0].keys())
        dict_writer.writeheader()
        dict_writer.writerows(output_results_csv_list)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    collate_results()