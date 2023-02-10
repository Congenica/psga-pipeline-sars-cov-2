from os.path import join as join_path  # used to join FS paths and S3 URIs
from typing import List

from scripts.synthetic.generate_pipeline_results_files import get_expected_output_files_per_sample


# pylint: disable=unused-argument
def get_expected_output_files(output_path: str, sample_ids: List[str], sequencing_technology: str) -> List[str]:
    """
    Return a list of ALL expected output paths
    """
    output_files_per_sample = get_expected_output_files_per_sample(
        output_path=output_path,
        sample_ids_result_files=sample_ids,
    )
    # generate a unified list of paths as non-sample results files must also be included
    output_files = [f["file"] for sample_files in output_files_per_sample.values() for f in sample_files]

    output_files.extend([join_path(output_path, "logs", f) for f in ["check_metadata.log"]])

    output_files.append(join_path(output_path, "results.csv"))
    output_files.append(join_path(output_path, "resultfiles.json"))

    return output_files
