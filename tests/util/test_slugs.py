import pytest

from app.scripts.util.slugs import get_file_with_type, FileType


@pytest.mark.parametrize("output_path", ["/", "s3://my_bucket/"])
@pytest.mark.parametrize("inner_dirs", [["a"], ["a", "b"], ["a", "b", "c"]])
@pytest.mark.parametrize(
    "filetypes",
    [
        [FileType(".bam", "bam/final"), FileType(".bam.bai", "bai/final")],
        [FileType(".fastq", "fastq/cleaned")],
    ],
)
@pytest.mark.parametrize("sample_id", ["roger", "clio"])
@pytest.mark.jira(identifier="aa8c8f26-fd37-41ef-b003-7aec00b71728", confirms="PSG-3621")
def test_get_file_with_type(
    output_path: str,
    inner_dirs: list[str],
    filetypes: list[FileType],
    sample_id: str,
):
    sample_files = get_file_with_type(
        output_path=output_path,
        inner_dirs=inner_dirs,
        filetypes=filetypes,
        sample_id=sample_id,
    )

    expected_files = []
    for ft in filetypes:
        expected_files.append(
            {
                "file": f"{output_path}{'/'.join(inner_dirs)}/{sample_id}{ft.extension}",
                "type": f"{ft.filetype}",
            }
        )
    assert sample_files == expected_files
