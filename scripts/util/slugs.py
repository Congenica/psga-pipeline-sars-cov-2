from os.path import join as join_path  # used to join FS paths and S3 URIs
from dataclasses import dataclass, field
from typing import Optional

SAMPLE_FILE_TYPE = dict[str, str]
RESULTFILES_TYPE = dict[str, list[SAMPLE_FILE_TYPE]]


@dataclass
class FileType:
    extension: str = field(metadata={"required": True})
    filetype: str = field(metadata={"required": True})
    order: Optional[int] = None


def get_file_with_type(
    output_path: str,
    inner_dirs: list[str],
    filetypes: list[FileType],
    sample_id: str,
) -> list[SAMPLE_FILE_TYPE]:
    """
    Return a dictionary { "file": "path/to/sample_id.ext", "type": "filetype", "order": None }
    """
    sample_result_files = []
    for result_file in filetypes:
        sample_result_file = {
            "file": join_path(output_path, *inner_dirs, f"{sample_id}{result_file.extension}"),
            "type": f"{result_file.filetype}",
        }
        if result_file.order:
            sample_result_file["order"] = result_file.order
        sample_result_files.append(sample_result_file)
    return sample_result_files
