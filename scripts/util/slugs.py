from os.path import join as join_path  # used to join FS paths and S3 URIs
from typing import Dict, List
from dataclasses import dataclass, field

SAMPLE_FILE_TYPE = Dict[str, str]
RESULTFILES_TYPE = Dict[str, List[SAMPLE_FILE_TYPE]]


@dataclass
class FileType:
    extension: str = field(metadata={"required": True})
    filetype: str = field(metadata={"required": True})


def get_file_with_type(
    output_path: str,
    inner_dirs: List[str],
    filetypes: List[FileType],
    sample_id: str,
) -> List[SAMPLE_FILE_TYPE]:
    """
    Return a dictionary { "file": "path/to/sample_id.ext", "type": "filetype" }
    """
    sample_files = [
        {
            "file": join_path(output_path, *inner_dirs, f"{sample_id}{elem.extension}"),
            "type": f"{elem.filetype}",
        }
        for elem in filetypes
    ]
    return sample_files
